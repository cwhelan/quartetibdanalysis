/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2011 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.quartetibdanalysis.tools;


import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.Sample;
import org.broadinstitute.hellbender.utils.samples.SampleDB;
import org.broadinstitute.hellbender.utils.samples.SampleDBBuilder;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

/**
 * Determines Identical-By-Descent regions in two siblings
 *
 * IBD0 regions - Siblings share neither parental haplotype
 * IBD1 regions - Siblings share one parental haplotype
 * IBD2 regions - Siblings share two parental haplotypes
 *
 * This tool implements two different methods for IBD region detection, which
 * should be used depending on whether the parental genotypes are available
 * (ie in a quartet sequencing project). Users can specify which type of analysis
 * to perform using the useParentalGenotypes parameter. If set to true, the tool
 * uses an HMM to label IBD regions based on the informativeness of each site (see
 * IDBObservation for details). If useParentalGenotypes is set to false, the tool
 * looks for sites at which minor alleles are shared between siblings and,
 * based on the level of sharing in sliding windows along the genome executes
 * a k-means model to determine the most likely IBD state combined
 * with median filters for smoothing. The latter method (based only on rare allele sharing)
 * is experimental and not recommended for general use.
 *
 * The tool scans the input pedigree file to detect sibling pairs, and computes
 * IBD regions for each pair of siblings detected in the pedigree file. We recommend
 * processing only one or two families at a time.
 *
 * IBD states are written to the output text file with the following columns:
 *
 * Sibling Pair ID: [Sibling Sample ID #1]-[Sibling Sample ID #2]
 * Chromosome
 * Region Start
 * Region End
 * IBD State: One of ZERO, ONE, or TWO
 *
 */
@CommandLineProgramProperties(
        summary = "Compute IBD blocks between siblings",
        oneLineSummary = "Compute IBD blocks between siblings",
        programGroup = QuartetIBDAnalysisProgramGroup.class
)
public class SiblingIBD extends VariantWalker {

    @Argument(fullName="mask", shortName = "mask", doc="BED file of regions to mask variants in")
    public FeatureInput<BEDFeature> mask = null;

    @Argument(shortName = "doNotUseParentalGenotypes",optional = false,fullName = "doNotUseParentalGenotypes", doc="Don't use parental genotypes, just count minor allele sharing between two siblings")
    private boolean doNotUseParentalGenotypes = false;

    @Argument(shortName = "gqThreshold",optional = true,fullName = "genotypeQualityThreshold", doc="Threshold on minimum GQ to include a site in the IBD calculation")
    private Integer gqThreshold = 50;

    @Argument(shortName = "hetABThreshold",optional = true,fullName = "hetABThreshold", doc="Threshold on minimum AB to include a het site in the IBD calculation")
    private Double hetABThreshold = 0.2;

    @Argument(shortName = "slidingWindowSize",optional = true,fullName = "slidingWindowSize", doc="Size of the sliding window to examine for shared alleles (Rare Allele Sharing Method Only)")
    private Integer slidingWindowSize = 1000;

    @Argument(shortName = "medianFilterSize",optional = true,fullName = "medianFilterSize", doc="Size of median filter to use to clean up IBD class predictions (Rare Allele Sharing Method Only)")
    private Integer medianFilterSize = 500;

    @Argument(shortName = "ibdRegionsFile", optional = true, doc = "File to which IBD regions should be written")
    private PrintStream ibdRegionsFile = null;

    @Argument(shortName = "useOriginalAF",optional = true,fullName = "useOriginalAF", doc="Use the AF_Orig field rather than the AF field for allele frequency (Rare Allele Sharing Method Only)")
    private boolean useOriginalAF = false;

    @Argument(shortName = "spFile",optional = true,fullName = "spFile", doc="File to output IBD observations at every unfiltered site for debugging purposes (Quartet Method Only)")
    private PrintStream spFile = null;

    @Argument(shortName = "countsFile",optional = true,fullName = "countsFile", doc="File to output sliding window counts for debugging purposes (Rare Allele Sharing Method Only)")
    private PrintStream countsFile = null;

    @Argument(shortName = "unfilteredIbdClassFile",optional = true,fullName = "unfilteredIbdClassFile", doc="File to output unfiltered IBD classes for debugging purposes (Rare Allele Sharing Method Only)")
    private PrintStream unfilteredIbdClassFile = null;

    @Argument(shortName = "filteredIbdClassFile",optional = true,fullName = "filteredIbdClassFile", doc="File to output filtered IBD classes for debugging purposes (Rare Allele Sharing Method Only)")
    private PrintStream filteredIbdClassFile = null;

    @Argument(shortName = "cnvCallsFile",optional = true,fullName = "cnvCallsFile", doc="VCF file containing CNV calls")
    private FeatureInput<VariantContext> cnvCalls = null;

    @Argument(fullName="pedigree", shortName="ped", doc="Pedigree file", optional=false)
    private File pedigreeFile = null;


    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written")
    public File outFile = null;

    private List<SiblingPair> siblingPairs;
    private Set<String> sampleNames;

    private IBDStateModel model = null;

    private VariantContextWriter vcfWriter = null;
    private SampleDB sampleDB;
    private Map<SiblingPair, IBDRegionAccumulator> ibdRegionAccumulators;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        final Map<String, VCFHeader> vcfHeaders = Collections.singletonMap(getDrivingVariantsFeatureInput().getName(), getHeaderForVariants());
        final SortedSet<String> vcfSamples = VcfUtils.getSortedSampleSet(vcfHeaders, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        sampleDB = initializeSampleDB();
        final Map<String, Set<Sample>> families = sampleDB.getFamilies(vcfSamples);
        sampleNames = new HashSet<>();
        final Collection<Set<Sample>> familyMemberSets = families.values();
        for(Set<Sample> familyMemberSet : familyMemberSets) {
            for (Sample sample : familyMemberSet) {
                sampleNames.add(sample.getID());
            }
        }

        siblingPairs = initSiblingPairs(vcfSamples);

        if (doNotUseParentalGenotypes) {
            model = new SiblingRareAlleleSharingModel(siblingPairs, gqThreshold, countsFile, unfilteredIbdClassFile, filteredIbdClassFile, useOriginalAF, slidingWindowSize, medianFilterSize);
        } else {
            model = new QuartetIBDStateModel(siblingPairs, gqThreshold, hetABThreshold, spFile);
        }

        final Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfHeaders.values(), true);

        headerLines.add(new VCFFormatHeaderLine("IBDS", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "IBD Max Likelihood State"));
        headerLines.add(new VCFFormatHeaderLine("IBDQ", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "IBD State Posteriors"));
        headerLines.add(new VCFFormatHeaderLine("IBDOBS", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "Potential Matching IBD States at SNP"));
        headerLines.add(new VCFFormatHeaderLine("IBDPA", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "IBD Parental Agreement"));
        headerLines.add(new VCFHeaderLine("source", "SiblingIBD"));

        if (outFile != null) {
            final SAMSequenceDictionary referenceDictionary = getReferenceDictionary();
            vcfWriter = createVCFWriter(referenceDictionary);
            vcfWriter.writeHeader(new VCFHeader(headerLines, sampleNames));
        }

        if (ibdRegionsFile != null) {
            ibdRegionsFile.print("PAIR_NAME" + "\t" + "SIB1" + "\t" + "SIB2" + "\t" + "CHR" + "\t" + "START" + "\t" + "END" + "\t" + "STATE" + "\t" + "PATERNAL_IND" + "\t" + "MATERNAL_IND" + "\t" + "SITES" + "\t" + "SITE_AGREEMENTS" + "\n");
        }

        ibdRegionAccumulators = new HashMap<>();
        for (final SiblingPair siblingPair : siblingPairs) {
            ibdRegionAccumulators.put(siblingPair, new IBDRegionAccumulator());
        }

    }

    private VariantContextWriter createVCFWriter(final SAMSequenceDictionary referenceDictionary) {
        VariantContextWriterBuilder vcWriterBuilder = new VariantContextWriterBuilder()
                                                            .clearOptions()
                                                            .setOutputFile(outFile);

        if (null != referenceDictionary) {
            vcWriterBuilder = vcWriterBuilder.setReferenceDictionary(referenceDictionary);
        }
        if (true) {
            vcWriterBuilder = vcWriterBuilder.setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        }
        for (Options opt : new Options[]{}) {
            vcWriterBuilder = vcWriterBuilder.setOption(opt);
        }

        return vcWriterBuilder.build();
    }

    private SampleDB initializeSampleDB() {
        final SampleDBBuilder sampleDBBuilder = new SampleDBBuilder(PedigreeValidationType.STRICT);
        sampleDBBuilder.addSamplesFromPedigreeFiles(Collections.singletonList(pedigreeFile));
        return sampleDBBuilder.getFinalSampleDB();
    }

    private List<SiblingPair> initSiblingPairs(final Set<String> vcfSamples) {
        System.out.println("init sibling pairs");
        final List<SiblingPair> siblingPairs = new ArrayList<>();
        final Map<String,Set<Sample>> families = sampleDB.getFamilies(vcfSamples);
        for(final String familyName : families.keySet()){
            System.out.println("Family: " + familyName);
            final Set<Sample> family = families.get(familyName);
            final Map<String, List<Sample>> familyMembersByParents = new HashMap<>();

            // first build up the list of family members by parents
            addFamilyMembersByParentsForFamily(vcfSamples, family, familyMembersByParents);

            // now find all sibling pairs in the family
            addSiblingPairsFromFamily(siblingPairs, familyMembersByParents);
        }

        return siblingPairs;
    }

    private void addSiblingPairsFromFamily(final List<SiblingPair> siblingPairs, final Map<String, List<Sample>> familyMembersByParents) {
        for (final String parentPairId : familyMembersByParents.keySet()) {
            System.out.println("parent pair: "+ parentPairId);
            final List<Sample> siblings = familyMembersByParents.get(parentPairId);
            System.out.println("siblings: " + siblings.size());
            if (siblings.size() >= 2) {
                for (int i  = 0; i < siblings.size(); i++) {
                    for (int j = i + 1; j < siblings.size(); j++) {
                        final Sample sib1 = siblings.get(i).getID().compareTo(siblings.get(j).getID()) < 0 ? siblings.get(i) : siblings.get(j);
                        final Sample sib2 = siblings.get(i).getID().compareTo(siblings.get(j).getID()) < 0 ? siblings.get(j) : siblings.get(i);
                        System.out.println("sib1: " + sib1.getID());
                        System.out.println("sib2: " + sib2.getID());
                        final SiblingPair siblingPair = new SiblingPair(sib1, sib2);
                        siblingPairs.add(siblingPair);
                    }
                }
            }
        }
    }

    private void addFamilyMembersByParentsForFamily(final Set<String> vcfSamples, final Set<Sample> family, final Map<String, List<Sample>> familyMembersByParents) {
        for(final Sample familyMember : family) {
            System.out.println("Processing family member " + familyMember);
            if (vcfSamples.contains(familyMember.getID())) {
                System.out.println("In vcf samples");
                final List<Sample> parents = sampleDB.getParents(familyMember);

                if (parents != null) {
                    System.out.println("Parents not null");
                    int validParents = 0;
                    for (Sample sample : parents) {
                        System.out.println("Parent: "+ sample);
                        if (!".".equals(sample.getID())) {
                            validParents += 1;
                        }
                    }
                    System.out.println("Valid parents: "+ validParents);
                    if (validParents == 2) {
                        final String parentPairId = parents.get(0).getID() + "-" + parents.get(1).getID();
                        if (!familyMembersByParents.containsKey(parentPairId)) {
                            familyMembersByParents.put(parentPairId, new ArrayList<>());
                        }
                        familyMembersByParents.get(parentPairId).add(familyMember);
                    }
                }
            }
        }
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if (!vc.isBiallelic()) {
            return;
        }

        if (!vc.isSNP()) {
            return;
        }

        if (vc.isFiltered()) {
            return;
        }

        if (featureContext.getValues(mask) != null && featureContext.getValues(mask).size() > 0) {
            return;
        }

        if (doNotUseParentalGenotypes) {
            final double af = SiblingRareAlleleSharingModel.getAF(useOriginalAF, vc);

            if (af > .4 && af < .6) {
                return;
            }
        }

        final VariantContext subsetVariantContext = vc.subContextFromSamples(sampleNames, false);

        final GenotypesContext genotypes = subsetVariantContext.getGenotypes(sampleNames);
        boolean hasAltGenotype = false;
        for (final String sample : sampleNames) {
            if (! genotypes.get(sample).isHomRef()) {
                hasAltGenotype = true;
                break;
            }
        }
        if (! hasAltGenotype) {
            return;
        }

        final List<VariantContext> overlappingCNVs = featureContext.getValues(cnvCalls);

        final List<VariantIBDState> result = new ArrayList<>();
        final SortedMap<Integer, List<VariantIBDState>> ibdStateModifications = new TreeMap<>();
        for (final SiblingPair siblingPair : siblingPairs) {
            final Genotype sib1Gt = genotypes.get(siblingPair.sib1.getID());
            final Genotype sib2Gt = genotypes.get(siblingPair.sib2.getID());

            if (sib1Gt.isCalled() && sib2Gt.isCalled()) {
                final List<VariantIBDState> variantIDBStatesForSibPair = model.addSite(subsetVariantContext,
                        genotypes, siblingPair, sib1Gt, sib2Gt, overlappingCNVs,outFile != null);
                if (variantIDBStatesForSibPair.size() == 0) {
                    continue;
                }
                for (final VariantIBDState variantIBDState : variantIDBStatesForSibPair) {
                    if (! ibdStateModifications.containsKey(variantIBDState.pos)) {
                        ibdStateModifications.put(variantIBDState.pos, new ArrayList<>());
                    }
                    ibdStateModifications.get(variantIBDState.pos).add(variantIBDState);
                }
                result.addAll(variantIDBStatesForSibPair);
            }
        }

        if (ibdStateModifications.keySet().size() == 0) {
            return;
        }

        for (final int variantIDBLoc : ibdStateModifications.keySet()) {

            final List<VariantIBDState> variantIBDStates = ibdStateModifications.get(variantIDBLoc);
            if (variantIBDStates.size() == 0) {
                continue;
            }
            if (outFile != null) {
                final VariantContext modifiedVC = variantIBDStates.get(0).vc;
                final VariantContextBuilder builder = new VariantContextBuilder(modifiedVC);
                final VariantContext newVC = addIBDStatesToVC(builder, modifiedVC, variantIBDStates);
                vcfWriter.add(newVC);
            }
        }

        if (logger.isDebugEnabled() && !result.isEmpty()) {
            logger.debug("emitting " + result);
        }
        reduce(result, ibdRegionAccumulators);

    }

    private VariantContext addIBDStatesToVC(final VariantContextBuilder builder, final VariantContext vc, final List<VariantIBDState> variantIBDStates) {
        final Map<String, String> ibdStateAttributes = new HashMap<>();
        final Map<String, String> ibdPosteriorAttributes = new HashMap<>();
        final Map<String, String> ibdObservationAttributes = new HashMap<>();
        final Map<String, String> ibdParentalAgreementAttributes = new HashMap<>();

        for (final VariantIBDState variantIBDState : variantIBDStates) {
            final SiblingPair siblingPair = variantIBDState.siblingPair;
            final Sample s1 = siblingPair.sib1;
            final Sample s2 = siblingPair.sib2;

            addIBDAttributesToMaps(s1, ibdStateAttributes, ibdPosteriorAttributes, ibdObservationAttributes, ibdParentalAgreementAttributes, variantIBDState);
            addIBDAttributesToMaps(s2, ibdStateAttributes, ibdPosteriorAttributes, ibdObservationAttributes, ibdParentalAgreementAttributes, variantIBDState);
        }

        final GenotypesContext newContext = GenotypesContext.create();
        for (final Genotype gt : vc.getGenotypes()) {
            GenotypeBuilder gtBuilder = new GenotypeBuilder(gt);
            if (ibdStateAttributes.containsKey(gt.getSampleName())) {
                gtBuilder = gtBuilder.attribute("IBDS", ibdStateAttributes.get(gt.getSampleName()));
            }
            if (ibdPosteriorAttributes.containsKey(gt.getSampleName())) {
                gtBuilder = gtBuilder.attribute("IBDQ", ibdPosteriorAttributes.get(gt.getSampleName()));
            }
            if (ibdObservationAttributes.containsKey(gt.getSampleName())) {
                gtBuilder = gtBuilder.attribute("IBDOBS", ibdObservationAttributes.get(gt.getSampleName()));
            }
            if (ibdParentalAgreementAttributes.containsKey(gt.getSampleName())) {
                gtBuilder = gtBuilder.attribute("IBDPA", ibdParentalAgreementAttributes.get(gt.getSampleName()));
            }

            newContext.add(gtBuilder.make());
        }
        return builder.genotypes(newContext).make();
    }

    private void addIBDAttributesToMaps(final Sample sample,
                                        final Map<String, String> ibdStateAttributes,
                                        final Map<String, String> ibdPosteriorAttributes,
                                        final Map<String, String> ibdObservationAttributes,
                                        final Map<String, String> ibdPaternalAgreementAttributes,
                                        final VariantIBDState variantIBDState) {
        final String ibdPairStateString = variantIBDState.siblingPair.sib1.getID() + "-" + variantIBDState.siblingPair.sib2.getID() + "=" + variantIBDState.ibdState.ordinal();
        final String sampleId = sample.getID();
        if (ibdStateAttributes.containsKey(sampleId)) {
            ibdStateAttributes.put(sampleId, ibdStateAttributes.get(sampleId) + ";" + ibdPairStateString);
        } else {
            ibdStateAttributes.put(sampleId, ibdPairStateString);
        }
        if (variantIBDState.statePosteriors != null) {
            final String ibdPosteriorString = variantIBDState.siblingPair.sib1.getID() + "-" + variantIBDState.siblingPair.sib2.getID() + "=" + formatPosteriors(variantIBDState);
            if (ibdPosteriorAttributes.containsKey(sampleId)) {
                ibdPosteriorAttributes.put(sampleId, ibdPosteriorAttributes.get(sampleId) + ";" + ibdPosteriorString  );
            } else {
                ibdPosteriorAttributes.put(sampleId, ibdPosteriorString);
            }
        }

        final String ibdPairObservationString = variantIBDState.siblingPair.sib1.getID() + "-" + variantIBDState.siblingPair.sib2.getID() + "=" + variantIBDState.observation.toString();
        if (ibdObservationAttributes.containsKey(sampleId)) {
            ibdObservationAttributes.put(sampleId, ibdObservationAttributes.get(sampleId) + ";" + ibdPairObservationString);
        } else {
            ibdObservationAttributes.put(sampleId, ibdPairObservationString);
        }

        final String ibdParentalAgreementString = variantIBDState.siblingPair.sib1.getID() + "-" + variantIBDState.siblingPair.sib2.getID() + "=" + (variantIBDState.parentalAgreement.toString());
        if (ibdPaternalAgreementAttributes.containsKey(sampleId)) {
            ibdPaternalAgreementAttributes.put(sampleId, ibdPaternalAgreementAttributes.get(sampleId) + ";" + ibdParentalAgreementString);
        } else {
            ibdPaternalAgreementAttributes.put(sampleId, ibdParentalAgreementString);
        }
    }

    private String formatPosteriors(final VariantIBDState variantIBDState) {
        final StringBuilder builder = new StringBuilder();
        for (int i = 0; i < variantIBDState.statePosteriors.length; i++) {
            final double p = variantIBDState.statePosteriors[i];
            builder.append(new Double(p));
            if (i < variantIBDState.statePosteriors.length - 1) {
                builder.append(",");
            }
        }
        return builder.toString();
    }


    public Map<SiblingPair, IBDRegionAccumulator> reduce(final List<VariantIBDState> valueList, final Map<SiblingPair, IBDRegionAccumulator> accumulators) {
        if (valueList != null) {
            accumulateValues(valueList, accumulators);
        }
        return accumulators;
    }

    private void accumulateValues(final List<VariantIBDState> valueList, final Map<SiblingPair, IBDRegionAccumulator> accumulators) {
        if (ibdRegionsFile != null) {
            for (final VariantIBDState value : valueList) {
                if (value != null) {
                    final SiblingPair siblingPair = value.siblingPair;
                    final IBDRegion region = accumulators.get(siblingPair).regionChange(value.chr, value.pos, value.ibdState, value.parentalAgreement, value.agrees);
                    if (region != null) {
                        if (region.state == IBDState.ONE) {
                            ibdRegionsFile.print(siblingPair.getName() + "\t" + siblingPair.sib1.getID() + "\t" + siblingPair.sib2.getID() + "\t" + region.chr + "\t" + region.start + "\t" + region.end + "\t" + region.state + "\t" + region.paternalIndicators + "\t" + region.maternalIndicators + "\t" + region.sites + "\t" + region.siteAgreements + "\n");
                        } else {
                            ibdRegionsFile.print(siblingPair.getName() + "\t" + siblingPair.sib1.getID() + "\t" + siblingPair.sib2.getID() + "\t" + region.chr + "\t" + region.start + "\t" + region.end + "\t" + region.state + "\t" + "NA" + "\t" + "NA" + "\t" + region.sites + "\t" + region.siteAgreements + "\n");
                        }
                    }
                }
            }
        }
    }


    @Override
    public Object onTraversalSuccess() {
        final SortedMap<Integer, List<VariantIBDState>> ibdStateModifications = new TreeMap<>();
        final List<VariantIBDState> finalValues = model.finalizeModel(vcfWriter != null);
        for (final VariantIBDState variantIBDState : finalValues) {
            if (! ibdStateModifications.containsKey(variantIBDState.pos)) {
                ibdStateModifications.put(variantIBDState.pos, new ArrayList<>());
            }
            ibdStateModifications.get(variantIBDState.pos).add(variantIBDState);
        }

        for (final int variantIDBLoc : ibdStateModifications.keySet()) {

            final List<VariantIBDState> variantIBDStates = ibdStateModifications.get(variantIDBLoc);

            if (variantIBDStates.size() == 0) {
                continue;
            }
            if (vcfWriter != null) {
                final VariantContext modifiedVC = variantIBDStates.get(0).vc;
                final VariantContextBuilder builder = new VariantContextBuilder(modifiedVC);
                final VariantContext newVC = addIBDStatesToVC(builder, modifiedVC, variantIBDStates);
                vcfWriter.add(newVC);
            }
        }

        accumulateValues(finalValues, ibdRegionAccumulators);
        if (ibdRegionsFile != null) {
            for (final SiblingPair siblingPair : ibdRegionAccumulators.keySet()) {
                final IBDRegionAccumulator accumulator = ibdRegionAccumulators.get(siblingPair);
                if (accumulator.chr == null) {
                    continue;
                }
                final IBDRegion region = accumulator.getFinalRegion();

                if (region.state == IBDState.ONE) {
                    ibdRegionsFile.print(siblingPair.getName() + "\t" + siblingPair.sib1.getID() + "\t" + siblingPair.sib2.getID() + "\t" + region.chr + "\t" + region.start + "\t" + region.end + "\t" + region.state + "\t" + region.paternalIndicators + "\t" + region.maternalIndicators + "\t" + region.sites + "\t" + region.siteAgreements + "\n");
                } else {
                    ibdRegionsFile.print(siblingPair.getName() + "\t" + siblingPair.sib1.getID() + "\t" + siblingPair.sib2.getID() + "\t" + region.chr + "\t" + region.start + "\t" + region.end + "\t" + region.state + "\t" + "NA" + "\t" + "NA" + "\t" + region.sites + "\t" + region.siteAgreements + "\n");
                }
            }
        }
        if (vcfWriter != null) {
            vcfWriter.close();
        }
        return null;
    }


}
