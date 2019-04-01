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

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.log4j.Logger;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by cwhelan on 5/27/14.
 *
 * The model strategy for the rare-allele sharing IBD prediction method.
 *
 */
public class SiblingRareAlleleSharingModel implements IBDStateModel {

    final protected static Logger logger = Logger.getLogger(SiblingRareAlleleSharingModel.class);
    private final List<SiblingPair> siblingPairs;
    private final Integer slidingWindowSize;
    private final Integer medianFilterSize;

    private final int gqThreshold;
    private PrintStream countsFile = null;
    private PrintStream unfilteredIbdClassFile = null;
    private PrintStream filteredIbdClassFile = null;
    private final boolean useOriginalAF;

    private final Map<SiblingPair, SlidingWindow<SharedMinorAlleleClass>> slidingWindowMap = new HashMap<>();
    private final Map<SiblingPair, IBDMedianFilter> ibdClustFilters = new HashMap<>();
    private Map<Integer, VariantContext> originalVariantContexts = new HashMap<>();

    private String currentContig;

    public SiblingRareAlleleSharingModel(final List<SiblingPair> siblingPairs, final int gqThreshold, final PrintStream countsFile, final PrintStream unfilteredIbdClassFile, final PrintStream filteredIbdClassFile, final boolean useOriginalAF, final Integer slidingWindowSize, final Integer medianFilterSize) {
        this.gqThreshold = gqThreshold;
        this.countsFile = countsFile;
        this.unfilteredIbdClassFile = unfilteredIbdClassFile;
        this.filteredIbdClassFile = filteredIbdClassFile;
        this.useOriginalAF = useOriginalAF;
        this.siblingPairs = siblingPairs;
        this.slidingWindowSize = slidingWindowSize;
        this.medianFilterSize = medianFilterSize;

        initSlidingWindows(siblingPairs);

        if (countsFile != null) {
            countsFile.print("CHROM\tPOS\tPAIR");
            for (final SharedMinorAlleleClass c : SharedMinorAlleleClass.values()) {
                countsFile.print("\t" + c.name());
            }
            countsFile.print("\n");
        }

        if (unfilteredIbdClassFile != null) {
            unfilteredIbdClassFile.print("CHROM\tPOS\tPAIR\tCLASS\n");
        }

        if (filteredIbdClassFile != null) {
            filteredIbdClassFile.print("CHROM\tPOS\tPAIR\tCLASS\n");
        }

    }

    static double getAF(final boolean useOriginalAF, final VariantContext vc) {
        if (useOriginalAF) {
            return vc.getAttributeAsDouble("AF_Orig", -1);
        } else {
            return vc.getAttributeAsDouble("AF", -1);
        }
    }

    @Override
    public List<VariantIBDState> addSite(final VariantContext vc, final GenotypesContext genotypes, final SiblingPair siblingPair, final Genotype sib1Gt, final Genotype sib2Gt, final List<VariantContext> overlappingCNVs, final boolean keepVc) {
        final List<VariantIBDState> result = new ArrayList<>();

        if (currentContig == null) {
            currentContig = vc.getContig();
        } else if (!currentContig.equals(vc.getContig())) {
            currentContig = vc.getContig();
            initSlidingWindows(siblingPairs);
        }

        final double af = getAF(useOriginalAF, vc);
        final List<Allele> sib1MinorAlleles = getMinorAlleles(af, sib1Gt);
        final List<Allele> sib2MinorAlleles = getMinorAlleles(af, sib2Gt);

        final int max = Math.max(sib1MinorAlleles.size(), sib2MinorAlleles.size());
        for (final Iterator<Allele> iterator = sib1MinorAlleles.iterator(); iterator.hasNext(); ) {
            final Allele a1 = iterator.next();
            if (sib2MinorAlleles.contains(a1)) {
                sib2MinorAlleles.remove(a1);
                iterator.remove();
            }
        }

        final int matches = max - Math.max(sib1MinorAlleles.size(), sib2MinorAlleles.size());
        final SharedMinorAllele sharedMinorAllele = new SharedMinorAllele(matches, max);
        if (!sharedMinorAllele.sharedAlleleClass().equals(SharedMinorAlleleClass.HOMREF_HOMREF) && gqThreshold <= Math.min(sib1Gt.getGQ(), sib2Gt.getGQ())) {

            final SlidingWindow<SharedMinorAlleleClass>.WindowCenterCount windowCenterCount = slidingWindowMap.get(siblingPair).next(vc.getStart(), sharedMinorAllele.sharedAlleleClass());
            if (windowCenterCount != null) {
                if (countsFile != null) {
                    countsFile.print(vc.getContig() + "\t" + windowCenterCount.start + "\t" + siblingPair.getName());
                    for (final SharedMinorAlleleClass c : SharedMinorAlleleClass.values()) {
                        countsFile.print("\t" + windowCenterCount.getCount(c));
                    }
                    countsFile.print("\n");
                }

                final AlleleCountKMeansModel model = new AlleleCountKMeansModel();
                final IBDState ibdClass = model.ibdClass(windowCenterCount);

                if (unfilteredIbdClassFile != null) {
                    unfilteredIbdClassFile.print(vc.getContig() + "\t" + windowCenterCount.start + "\t" + siblingPair.getName() + "\t" + ibdClass + "\n");
                }
                final IBDMedianFilter.FilteredValue filteredValue = ibdClustFilters.get(siblingPair).next(windowCenterCount.start, ibdClass);
                if (filteredValue != null) {
                    final IBDState filteredIbdClass = filteredValue.ibdClass;
                    final VariantContext variantContext = originalVariantContexts.get(filteredValue.start);
                    final VariantIBDState variantIBDState = new VariantIBDState(variantContext.getContig(),
                            variantContext.getStart(),
                            keepVc ? variantContext : null,
                            siblingPair, filteredIbdClass,
                            null, IBDParentalAgreement.UNKNOWN, false, null);
                    if (logger.isDebugEnabled()) {
                        logger.debug("adding to result: " + variantIBDState);
                    }
                    result.add(variantIBDState);
                    if (filteredIbdClassFile != null) {
                        if (logger.isDebugEnabled()) {
                            logger.debug("writing to fic file: " + vc.getContig() + "\t" + filteredValue.start + "\t" + siblingPair.getName() + "\t" + filteredIbdClass + "\n");
                        }
                        filteredIbdClassFile.print(vc.getContig() + "\t" + filteredValue.start + "\t" + siblingPair.getName() + "\t" + filteredIbdClass + "\n");
                    }
                }
            }
        }

        return result;
    }

    @Override
    public List<VariantIBDState> finalizeModel(final boolean keepVc) {
        return new ArrayList<>();
    }

    private void initSlidingWindows(final List<SiblingPair> siblingPairs) {
        for (final SiblingPair siblingPair : siblingPairs) {
            slidingWindowMap.put(siblingPair, new SlidingWindow<SharedMinorAlleleClass>(slidingWindowSize));
            ibdClustFilters.put(siblingPair, new IBDMedianFilter(medianFilterSize));
        }
        originalVariantContexts = new HashMap<>();
    }

    private List<Allele> getMinorAlleles(final double af, final Genotype gt) {
        final List<Allele> sib1MinorAlleles = new ArrayList<>();
        for (final Allele allele : gt.getAlleles()) {
            if (allele.isCalled() && (af < 0.5 ? allele.isNonReference() : !allele.isNonReference())) {
                sib1MinorAlleles.add(allele);
            }
        }
        return sib1MinorAlleles;
    }

}
