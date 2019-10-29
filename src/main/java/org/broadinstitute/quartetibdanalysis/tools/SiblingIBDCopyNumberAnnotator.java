package org.broadinstitute.quartetibdanalysis.tools;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.linear.*;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.samples.PedReader;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.SampleDB;
import org.broadinstitute.hellbender.utils.samples.SampleDBBuilder;
import org.broadinstitute.sv.annotation.GenotypedVariantFactory;
import org.broadinstitute.sv.util.GenomeInterval;
import org.broadinstitute.sv.util.IntervalTreeMap;
import org.broadinstitute.sv.util.io.ErrorCheckingPrintWriter;
import org.broadinstitute.sv.util.ped.PedFileParser;
import org.broadinstitute.sv.util.ped.Pedigree;
import org.broadinstitute.sv.util.vc.GenotypeUtilities;
import org.broadinstitute.sv.util.vc.VariantContextUtilities;
import org.broadinstitute.sv.util.vcf.VCFReader;

import java.io.*;
import java.util.*;

/**
 * Created by cwhelan on 7/8/15.
 */
public class SiblingIBDCopyNumberAnnotator {
    private static final Logger mLogger = Logger.getLogger(SiblingIBDCopyNumberAnnotator.class);
    public static final String M2_MODE_DOUBLE_IBD = "doubleIBD";
    public static final String M2_MODE_LIKELIHOOD = "likelihood";

    private GenotypedVariantFactory mVariantFactory = null;
    private Map<SiblingPair, IntervalTreeMap<IBDState>> ibdStateMaps;
    private List<SiblingPair> siblingPairs;
    //private Map<String, Pedigree> mPedigrees;
    private String[] samples;
    private Map<SiblingPair, String> mSibPairToFamilyMap;
    private Float minCallRate = null;
    private boolean processErrorProbs = true;
    private String m2Mode;

    private Map<String, VariantContext> ends1ById = new HashMap<>();
    private ErrorCheckingPrintWriter endsConfidenceWriter;
    private int endConfidenceLikelihoodThreshold = 20;
    private PrintWriter reportWriter;
    private SampleDB finalSampleDB;

    public void initialize(final boolean filterGenotypes,
                           final Double genotypeQualityThreshold,
                           final File siblingIBDFile,
                           final Float minCallRate,
                           final Integer endConfidenceLikelihoodThreshold,
                           final Boolean processErrorProbs,
                           final String m2Mode,
                           final File ends1VCF,
                           final File confidenceReport,
                           final List<File> pedigreeFileList,
                           final boolean writeReport,
                           final File outputFile) {
        mVariantFactory = new GenotypedVariantFactory();
        mVariantFactory.setFilterGenotypes(filterGenotypes);
        mVariantFactory.setGenotypeQualityThreshold(genotypeQualityThreshold);

        this.minCallRate = minCallRate;
        if (endConfidenceLikelihoodThreshold != null) {
            this.endConfidenceLikelihoodThreshold = endConfidenceLikelihoodThreshold;
        }
        this.processErrorProbs = processErrorProbs;

        this.m2Mode = m2Mode == null ? M2_MODE_DOUBLE_IBD : m2Mode;

        if (this.m2Mode.equals(M2_MODE_LIKELIHOOD)) {
            try (VCFReader ends1Reader = new VCFReader(ends1VCF)) {
                for (VariantContext vc : ends1Reader) {
                    ends1ById.put(vc.getID(), vc);
                }
            }

            endsConfidenceWriter = new ErrorCheckingPrintWriter(confidenceReport);
        }

        ibdStateMaps = new HashMap<>();
        siblingPairs = new ArrayList<>();
        List<String> samplesList = new ArrayList<>();



        List<File> pedigreeFiles = pedigreeFileList;
        System.out.println("about to make ped file parser");
        try {
            //PedFileParser pedFileParser = new PedFileParser();

            final SampleDBBuilder sampleDBBuilder = new SampleDBBuilder(PedigreeValidationType.STRICT);
            sampleDBBuilder.addSamplesFromPedigreeFiles(pedigreeFiles);
            finalSampleDB = sampleDBBuilder.getFinalSampleDB();


        } catch (Throwable t) {
            t.printStackTrace();
        }

        mSibPairToFamilyMap = new HashMap<>();

        try {
            if (!siblingIBDFile.isFile()) {
                throw new RuntimeException("File not found: " + siblingIBDFile);
            }
            BufferedReader reader = new BufferedReader(new FileReader(siblingIBDFile));
            String[] header = reader.readLine().split("\t");

            while (true) {
                String line = reader.readLine();
                if (line == null) {
                    break;
                }
                String[] fields = line.split("\t");
                SiblingPair siblingPair = new SiblingPair( fields[1], fields[2]);

                Interval interval = new Interval(fields[3], Integer.parseInt(fields[4]), Integer.parseInt(fields[5]));
                String state = fields[6];
                String paternalIndicators = fields[7];
                String maternalIndicators = fields[8];
                String sites = fields[9];
                Integer sitesDouble = null;
                if (!"NA".equals(sites)) {
                    sitesDouble = Integer.parseInt(sites);
                }

                Double paternalIndicatorsDouble = null;
                if (!"NA".equals(paternalIndicators)) {
                    paternalIndicatorsDouble = Double.parseDouble(paternalIndicators);
                }
                Double maternalIndicatorsDouble = null;
                if (!"NA".equals(maternalIndicators)) {
                    maternalIndicatorsDouble = Double.parseDouble(maternalIndicators);
                }

                final IBDState ibdState = new IBDState(state, sitesDouble, paternalIndicatorsDouble, maternalIndicatorsDouble);

                if (! ibdStateMaps.containsKey(siblingPair)) {
                    ibdStateMaps.put(siblingPair, new IntervalTreeMap<IBDState>());
                    siblingPairs.add(siblingPair);
                    samplesList.add(siblingPair.sib1);
                    samplesList.add(siblingPair.sib2);
                }

                ibdStateMaps.get(siblingPair).put(interval, ibdState);
                final String familyId = findFamilyId(siblingPair);
                if (familyId == null) {
                    throw new GATKException("Couldn't find a family ID for sib pair " + siblingPair + " in the pedigree");
                }
                mSibPairToFamilyMap.put(siblingPair, familyId);

                if (! finalSampleDB.getFamilyIDs().contains(familyId)) {
                    throw new GATKException("Can't find familyId " + familyId + " in the pedigree");
                }
                samplesList.add(finalSampleDB.getSample(siblingPair.sib1).getPaternalID());
                samplesList.add(finalSampleDB.getSample(siblingPair.sib1).getMaternalID());
            }
        }  catch (IOException exc) {
            throw new RuntimeException(exc.getMessage(), exc);
        }
        samples = samplesList.toArray(new String[samplesList.size()]);

        if (writeReport) {
            reportWriter = getReportWriter(outputFile);
            writeReportHeader();
        }
    }

    private PrintWriter getReportWriter(final File outputFile) {
        try {
            return new ErrorCheckingPrintWriter(new FileWriter(outputFile));
        } catch (IOException exc) {
            throw new GATKException("Cannot open output file: " + outputFile + ": " + exc.getMessage(), exc);
        }
    }

    private String findFamilyId(final SiblingPair siblingPair) {
        final String sib1FamId = finalSampleDB.getSample(siblingPair.sib1).getFamilyID();
        final String sib2FamId = finalSampleDB.getSample(siblingPair.sib2).getFamilyID();

        if (sib1FamId.equals(sib2FamId)) {
            return sib1FamId;
        } else {
            return null;
        }
    }

    protected void addField(StringBuilder builder, Object value) {
        if (builder.length() > 0) {
            builder.append('\t');
        }
        builder.append((value == null) ? "NA" : value.toString());
    }

    protected String formatDouble(String format, Double value) {
        if (value == null || value.isNaN()) {
            return "NA";
        }
        return String.format(format, value);
    }

    // Centralize this to help standardize interval handling.
    protected GenomeInterval getVariantInterval(VariantContext vc) {
        return VariantContextUtilities.getInterval(vc);
    }


    private void writeReportHeader() {
        StringBuilder builder = new StringBuilder();
        addField(builder, "ID");
        addField(builder, "SIB_PAIR");
        addField(builder, "SIB1");
        addField(builder, "SIB2");
        addField(builder, "IBDSTATE");
        addField(builder, "SIB1_CN");
        addField(builder, "SIB1_CNQ");
        addField(builder, "SIB1_CNP");
        addField(builder, "SIB2_CN");
        addField(builder, "SIB2_CNQ");
        addField(builder, "SIB2_CNP");
        addField(builder, "PAT_ID");
        addField(builder, "MAT_ID");
        addField(builder, "PAT_CN");
        addField(builder, "PAT_CNQ");
        addField(builder, "PAT_CNP");
        addField(builder, "MAT_CN");
        addField(builder, "MAT_CNQ");
        addField(builder, "MAT_CNP");
        addField(builder, "IBDSTATE-REFINED");
        addField(builder, "PARENT_CONFIDENCE");
        addField(builder, "SITES");
        addField(builder, "DISTANCE_TO_IBD_STATE_CHANGE");
        addField(builder, "CONST_IBD0");
        addField(builder, "CONST_IBD1");
        addField(builder, "CONST_IBD1P");
        addField(builder, "CONST_IBD1M");
        addField(builder, "CONST_IBD2");
        addField(builder, "STATE_CONSISTENT");
        addField(builder, "STATE_CONSISTENT_MUTGAIN1_SIB1");
        addField(builder, "STATE_CONSISTENT_MUTGAIN1_SIB2");
        addField(builder, "STATE_CONSISTENT_MUTLOSS1_SIB1");
        addField(builder, "STATE_CONSISTENT_MUTLOSS1_SIB1");

        addField(builder, "TOTAL_PROB_IBD0");
        addField(builder, "TOTAL_PROB_IBD1");
        addField(builder, "TOTAL_PROB_IBD1P");
        addField(builder, "TOTAL_PROB_IBD1M");
        addField(builder, "TOTAL_PROB_IBD2");

        addField(builder, "ERROR_COND");
        reportWriter.println(builder.toString());
    }

    private void writeReportLine(VariantContext vc, SiblingPair siblingPair, IBDState ibdState, final Integer distanceToIBDStateChange, String errorCond) {
        final Genotype sib1Gt = vc.getGenotype(siblingPair.sib1);
        final Genotype sib2Gt = vc.getGenotype(siblingPair.sib2);
        final String family = mSibPairToFamilyMap.get(siblingPair);

        final String father = finalSampleDB.getSample(siblingPair.sib1).getPaternalID();
        final String mother = finalSampleDB.getSample(siblingPair.sib1).getMaternalID();
        final Genotype patGt = vc.getGenotype(father);
        final Genotype matGt = vc.getGenotype(mother);

        StringBuilder builder = new StringBuilder();
        addField(builder, vc.getID());
        addField(builder, siblingPair.getName());
        addField(builder, siblingPair.sib1);
        addField(builder, siblingPair.sib2);
        addField(builder, ibdState != null ? ibdState.state : "NA");
        final int sib1Cn = formatCNFields(builder, sib1Gt);
        final int sib2Cn = formatCNFields(builder, sib2Gt);
        addField(builder, father);
        addField(builder, mother);
        final int patCn = formatCNFields(builder, patGt);
        final int matCn = formatCNFields(builder, matGt);
        addField(builder, ibdState != null ? ibdState.getRefinedState() : "NA");
        addField(builder, ibdState != null ? ibdState.getParentConfidence() : "NA");

        addField(builder, ibdState != null ? ibdState.getSites() : "NA");

        addField(builder, distanceToIBDStateChange != null ? distanceToIBDStateChange : "NA");
        addField(builder, formatBoolean(consistentWithIDB0(sib1Cn, sib2Cn, patCn, matCn)));
        addField(builder, formatBoolean(consistentWithIDB1(sib1Cn, sib2Cn, patCn, matCn)));
        addField(builder, formatBoolean(consistentWithIDB1P(sib1Cn, sib2Cn, patCn, matCn)));
        addField(builder, formatBoolean(consistentWithIDB1M(sib1Cn, sib2Cn, patCn, matCn)));
        addField(builder, formatBoolean(consistentWithIDB2(sib1Cn, sib2Cn)));
        addField(builder, formatBoolean(stateConsistent(ibdState, sib1Cn, sib2Cn, patCn, matCn)));

        final int sib1GainCn = sib1Cn - 1;
        addField(builder, formatBoolean(stateConsistent(ibdState, sib1GainCn, sib2Cn, patCn, matCn)));

        final int sib2GainCn = sib2Cn - 1;
        addField(builder, formatBoolean(stateConsistent(ibdState, sib1Cn, sib2GainCn, patCn, matCn)));

        final int sib1LossCn = sib1Cn + 1;
        addField(builder, formatBoolean(stateConsistent(ibdState, sib1LossCn, sib2Cn, patCn, matCn)));

        final int sib2LossCn = sib2Cn + 1;
        addField(builder, formatBoolean(stateConsistent(ibdState, sib1Cn, sib2LossCn, patCn, matCn)));



        if (sib1Cn != -1 && sib2Cn != -1 && patCn != -1 && matCn != -1 && (errorCond == null || processErrorProbs)) {
            final double[] totalStateProbabilities = computeTotalStateProbabilities(sib1Gt, sib2Gt, patGt, matGt);

            addField(builder, formatDouble("%.4f", totalStateProbabilities[0]));
            addField(builder, formatDouble("%.4f", totalStateProbabilities[1]));
            addField(builder, formatDouble("%.4f", totalStateProbabilities[2]));
            addField(builder, formatDouble("%.4f", totalStateProbabilities[3]));
            addField(builder, formatDouble("%.4f", totalStateProbabilities[4]));
        } else {
            addField(builder, "NA");
            addField(builder, "NA");
            addField(builder, "NA");
            addField(builder, "NA");
            addField(builder, "NA");
        }


        addField(builder, errorCond);

        reportWriter.println(builder.toString());
    }

    private double[] computeTotalStateProbabilities(final Genotype sib1Gt, final Genotype sib2Gt, final Genotype patGt, final Genotype matGt) {
        // [0]: IBD0 [1]: IBD1 [2]: IBD1-P [3]: IBD1-M [4]: IBD2
        double[] totalProbs = new double[5];

        final double[] sib1GtPosteriors = GenotypeUtilities.getCopyNumberLogPosteriors(sib1Gt);
        final double[] sib2GtPosteriors = GenotypeUtilities.getCopyNumberLogPosteriors(sib2Gt);
        final double[] patGtPosteriors = GenotypeUtilities.getCopyNumberLogPosteriors(patGt);
        final double[] matGtPosteriors = GenotypeUtilities.getCopyNumberLogPosteriors(matGt);

        for (int i = 0; i < totalProbs.length; i++) {
            totalProbs[i] = Double.NEGATIVE_INFINITY;
        }

        for (int sib1Cn = 0; sib1Cn < sib1GtPosteriors.length; sib1Cn++) {
            for (int sib2Cn = 0; sib2Cn < sib2GtPosteriors.length; sib2Cn++) {
                for (int patCn = 0; patCn < patGtPosteriors.length; patCn++) {
                    for (int matCn = 0; matCn < matGtPosteriors.length; matCn++) {

                        if (consistentWithIDB0(sib1Cn, sib2Cn, matCn, patCn)) {
                            final double pGenotypeConfig = sib1GtPosteriors[sib1Cn] + sib2GtPosteriors[sib2Cn] + patGtPosteriors[patCn] + matGtPosteriors[matCn];
                            totalProbs[0] = MathUtils.log10sumLog10(new double[] {totalProbs[0], pGenotypeConfig});
                        }

                        if (consistentWithIDB1(sib1Cn, sib2Cn, matCn, patCn)) {
                            final double pGenotypeConfig = sib1GtPosteriors[sib1Cn] + sib2GtPosteriors[sib2Cn] + patGtPosteriors[patCn] + matGtPosteriors[matCn];
                            totalProbs[1] = MathUtils.log10sumLog10(new double[] {totalProbs[1], pGenotypeConfig});
                        }

                        if (consistentWithIDB1P(sib1Cn, sib2Cn, matCn, patCn)) {
                            final double pGenotypeConfig = sib1GtPosteriors[sib1Cn] + sib2GtPosteriors[sib2Cn] + patGtPosteriors[patCn] + matGtPosteriors[matCn];
                            totalProbs[2] = MathUtils.log10sumLog10(new double[] {totalProbs[2], pGenotypeConfig});
                        }

                        if (consistentWithIDB1M(sib1Cn, sib2Cn, matCn, patCn)) {
                            final double pGenotypeConfig = sib1GtPosteriors[sib1Cn] + sib2GtPosteriors[sib2Cn] + patGtPosteriors[patCn] + matGtPosteriors[matCn];
                            totalProbs[3] = MathUtils.log10sumLog10(new double[] {totalProbs[3], pGenotypeConfig});
                        }

                        if (consistentWithIDB2(sib1Cn, sib2Cn)) {
                            final double pGenotypeConfig = sib1GtPosteriors[sib1Cn] + sib2GtPosteriors[sib2Cn] + patGtPosteriors[patCn] + matGtPosteriors[matCn];
                            totalProbs[4] = MathUtils.log10sumLog10(new double[] {totalProbs[4], pGenotypeConfig});
                        }

                    }
                }
            }
        }
        return totalProbs;
    }


    private int formatCNFields(final StringBuilder builder, final Genotype gt) {
        if (gt == null) {
            addEmptyCNFields(builder);
            return -1;
        } else {
            final double[] copyNumberLogPosteriors = GenotypeUtilities.getCopyNumberLogPosteriors(gt);
            if (! GenotypeUtilities.isHardFiltered(gt)) {
                final Integer copyNumber = GenotypeUtilities.getCopyNumber(gt);
                final Double copyNumberQuality = GenotypeUtilities.getCopyNumberQuality(gt);
                if (copyNumber == null || copyNumberQuality == null) {
                    addEmptyCNFields(builder);
                    return -1;
                }
                addField(builder, copyNumber);
                addField(builder,  formatDouble("%.4f", copyNumberQuality));
                formatPosteriors(builder, copyNumberLogPosteriors);
                return copyNumber;
            } else {
                final Double copyNumberQuality = GenotypeUtilities.getCopyNumberQuality(gt);
                if (copyNumberQuality == null) {
                    addEmptyCNFields(builder);
                    return -1;
                }
                addField(builder, -1);
                addField(builder,  formatDouble("%.4f", copyNumberQuality));
                formatPosteriors(builder, copyNumberLogPosteriors);
                return -1;
            }
        }
    }

    private void addEmptyCNFields(final StringBuilder builder) {
        addField(builder, -1);
        addField(builder,  "NA");
        addField(builder,  "NA");
    }

    private void formatPosteriors(final StringBuilder builder, final double[] copyNumberLogPosteriors) {
        if (copyNumberLogPosteriors != null) {
            StringBuilder posteriorBuilder = new StringBuilder();
            for (int i = 0; i < copyNumberLogPosteriors.length; i++) {
                posteriorBuilder.append(formatDouble("%.4f", copyNumberLogPosteriors[i]));
                if (i < copyNumberLogPosteriors.length - 1) {
                    posteriorBuilder.append(",");
                }
            }
            addField(builder, posteriorBuilder.toString());
        } else {
            addField(builder, "NA");
        }
    }

    private boolean stateConsistent(final IBDState ibdState, final int sib1Gt, final int sib2Gt, final int patGt, final int matGt) {
        if (ibdState == null || ibdState.getRefinedState() == null) {
            return false;
        }
        final String refinedState = ibdState.getRefinedState();
        return ((refinedState.equals("ZERO") && consistentWithIDB0(sib1Gt, sib2Gt, patGt, matGt) != null && consistentWithIDB0(sib1Gt, sib2Gt, patGt, matGt)) ||
                        (refinedState.equals("ONE-P")  && consistentWithIDB1P(sib1Gt, sib2Gt, patGt, matGt) != null && consistentWithIDB1P(sib1Gt, sib2Gt, patGt, matGt)) ||
                        (refinedState.equals("ONE-M") && consistentWithIDB1M(sib1Gt, sib2Gt, patGt, matGt) != null && consistentWithIDB1M(sib1Gt, sib2Gt, patGt, matGt)) ||
                        (refinedState.equals("TWO") && consistentWithIDB2(sib1Gt, sib2Gt) != null && consistentWithIDB2(sib1Gt, sib2Gt)));
    }

    private String formatBoolean(final Boolean aBoolean) {
        return aBoolean == null ? "NA" : (aBoolean ? "TRUE" : "FALSE");
    }

    public int getSampleIdx(String sampleName, String[] samples) {
        for (int i = 0; i < samples.length; i++) {
            if (sampleName.equals(samples[i])) {
                return i;
            }
        }
        return -1;
    }

    public Boolean consistentWithIDB2(int sib1Gt, int sib2Gt) {
        if (sib1Gt == -1 || sib2Gt == -1) return null;
        return sib1Gt == sib2Gt;
    }

    public Boolean consistentWithIDB0(int sib1Gt, int sib2Gt, int patGt, int matGt) {
        if (sib1Gt == -1 || sib2Gt == -1 || patGt == -1 || matGt == -1) return null;
        return sib1Gt + sib2Gt == patGt + matGt;
    }

    public Boolean consistentWithIDB1(int sib1Gt, int sib2Gt, int patGt, int matGt) {
        if (sib1Gt == -1 || sib2Gt == -1 || patGt == -1 || matGt == -1) return null;
        if (patGt % 2 == 0 && matGt % 2 == 0 && (sib1Gt + sib2Gt) % 2 == 1) return false;
        if (patGt % 2 == 1 && matGt % 2 == 1 && (sib1Gt + sib2Gt) % 2 == 0) return false;
        return true;
    }

    public Boolean consistentWithIDB1P(int sib1Gt, int sib2Gt, int patGt, int matGt) {
        if (sib1Gt == -1 || sib2Gt == -1 || patGt == -1 || matGt == -1) return null;
        RealMatrix coefficients =
                new Array2DRowRealMatrix(new double[][] { { 1,1,0,0 }, { 0,0,1,1 }, { 1,0,1,0 }, {1,0,0,1} },
                        false);
        DecompositionSolver solver = new LUDecomposition(coefficients).getSolver();
        //double[] solution = solver.solve(new double[] {patGt,matGt,sib1Gt,sib2Gt});
        RealVector solution = solver.solve(new ArrayRealVector(new double[] {patGt,matGt,sib1Gt,sib2Gt} ));

        return isAllIntegers(solution.toArray());
    }

    public Boolean consistentWithIDB1M(int sib1Gt, int sib2Gt, int patGt, int matGt) {
        // pretend the mother is the father and this should work out
        return consistentWithIDB1P(sib1Gt, sib2Gt, matGt, patGt);
    }


    private Boolean isAllIntegers(final double[] solution) {
        for (int i = 0; i < solution.length; i++) {
            if (! (Math.abs(solution[i] - Math.rint(solution[i])) < 0.001)) {
                return false;
            }
        }
        return true;
    }


    public Map<String, Object> annotate(final VariantContext vc) {
        if (vc == null) {
            return null;
        }

        if (!vc.hasGenotypes()) {
            return null;
        }

        if (vc.isFiltered()) {
            // System.out.println("Variant is filtered, skipping");
            return null;
        }

        if (minCallRate != null && vc.getAttributeAsDouble("GSCALLRATE", 0.0) < minCallRate.doubleValue()) {
            return null;
        }

        //GenotypedVariant variant = mVariantFactory.createVariant(vc, samples);
        //VariantGenotypeModel gtModel = variant.getGenotypeModel();

        for (SiblingPair siblingPair : ibdStateMaps.keySet()) {
            GenomeInterval gi = getVariantInterval(vc);
            Interval interval;

            if (! isM2Genotyping(vc)) {
                interval = new Interval(gi.getSequenceName(), gi.getStart(), gi.getEnd());
                final Collection<IBDState> overlappingIBDRegions = ibdStateMaps.get(siblingPair).getOverlapping(interval);

                processSingleInterval(vc, siblingPair, gi, interval, overlappingIBDRegions, vc.getID(), null, null);
                continue;

            } else {
                if (m2Mode.equals(M2_MODE_DOUBLE_IBD)) {
                    processDoubleIntervals(vc, siblingPair, gi);
                    continue;

                } else if (m2Mode.equals(M2_MODE_LIKELIHOOD)) {
                    final String vcID = vc.getID();

                    final String[] idFields = vcID.split("_");
                    final String end1Id = "CNV" + "_" + idFields[3] + "_" + (Integer.parseInt(idFields[4]) - 1) + "_" + idFields[5];
                    final String end2Id = "CNV" + "_" + idFields[6] + "_" + (Integer.parseInt(idFields[7]) - 1) + "_" + idFields[8];

                    final VariantContext end1Vc = ends1ById.get(end1Id);
                    final VariantContext end2Vc = ends1ById.get(end2Id);

                    if (end1Vc == null) {
                        mLogger.warn("could not find end genotypes for " + end1Id);
                        return null;
                    }

                    if (end2Vc == null) {
                        mLogger.warn("could not find end genotypes for " + end2Id);
                        return null;
                    }

                    final String sib1 = siblingPair.sib1;
                    final String sib2 = siblingPair.sib2;

                    final String father = finalSampleDB.getSample(siblingPair.sib1).getPaternalID();
                    final String mother = finalSampleDB.getSample(siblingPair.sib1).getMaternalID();

                    if (! vc.hasGenotype(sib1) ||
                            ! vc.hasGenotype(sib2) ||
                            ! vc.hasGenotype(father) ||
                            ! vc.hasGenotype(mother)) {
                        mLogger.warn("incomplete family genotypes for " + vcID + "/" + siblingPair.getName());
                        return null;
                    }

                    if (! vc.getGenotype(sib1).hasExtendedAttribute("CNL") ||
                            ! vc.getGenotype(sib2).hasExtendedAttribute("CNL") ||
                            ! vc.getGenotype(father).hasExtendedAttribute("CNL") ||
                            ! vc.getGenotype(mother).hasExtendedAttribute("CNL")) {
                        mLogger.warn("missing genotype likelihoods for " + vcID + "/" + siblingPair.getName());
                        return null;
                    }

                    LikelihoodCombination bestSib1Likelihood = getBestLikelihood(vc, end1Vc, end2Vc, sib1);
                    LikelihoodCombination bestSib2Likelihood = getBestLikelihood(vc, end1Vc, end2Vc, sib2);

                    LikelihoodCombination bestFatherLikelihood = getBestLikelihood(vc, end1Vc, end2Vc, father);
                    LikelihoodCombination bestMotherLikelihood = getBestLikelihood(vc, end1Vc, end2Vc, mother);

                    if (bestSib1Likelihood == null || bestSib2Likelihood == null || bestFatherLikelihood == null || bestMotherLikelihood == null) {
                        mLogger.warn("null likelihoods for " + vcID + "/" + vcID);
                        return null;
                    }

                    final boolean allATwo = bestSib1Likelihood.end1Cn == 2 &&
                            bestSib2Likelihood.end1Cn == 2 &&
                            bestFatherLikelihood.end1Cn == 2 &&
                            bestMotherLikelihood.end1Cn == 2;
                    final boolean allBTwo = bestSib1Likelihood.end2Cn == 2 &&
                            bestSib2Likelihood.end2Cn == 2 &&
                            bestFatherLikelihood.end2Cn == 2 &&
                            bestMotherLikelihood.end2Cn == 2;

                    endsConfidenceWriter.println(
                            vcID + "\t" + siblingPair.getName() + "\t" +
                                    sib1 + "\t" + bestSib1Likelihood.end1Cn + "\t" + bestSib1Likelihood.end2Cn + "\t" + bestSib1Likelihood.sdCn + "\t" + bestSib1Likelihood.ll + "\t" +
                                    sib2 + "\t" + bestSib2Likelihood.end1Cn + "\t" + bestSib2Likelihood.end2Cn + "\t" + bestSib2Likelihood.sdCn + "\t" + bestSib2Likelihood.ll + "\t" +
                                    father + "\t" + bestFatherLikelihood.end1Cn + "\t" + bestFatherLikelihood.end2Cn + "\t" + bestFatherLikelihood.sdCn + "\t" + bestFatherLikelihood.ll + "\t" +
                                    mother + "\t" + bestMotherLikelihood.end1Cn + "\t" + bestMotherLikelihood.end2Cn + "\t" + bestMotherLikelihood.sdCn + "\t" + bestMotherLikelihood.ll);

                    String paralogId;
                    if (allATwo && ! allBTwo) {
                        paralogId = end2Id;
                    } else if (! allATwo && allBTwo) {
                        paralogId = end1Id;
                    } else {
                        mLogger.warn("Not a mix of all2 and notAll2 for " + vcID + "/" + siblingPair.getName());
                        return null;
                    }


                    final Double minSibLL = Math.min(bestSib1Likelihood.ll, bestSib2Likelihood.ll);
                    final Double minParentLL = Math.min(bestSib1Likelihood.ll, bestSib2Likelihood.ll);
                    final Double minFamilyLL = Math.min(minSibLL, minParentLL);

                    final String[] paralogIdFields = paralogId.split("_");
                    interval = new Interval(paralogIdFields[1], Integer.parseInt(paralogIdFields[2]),  Integer.parseInt(paralogIdFields[3]));
                    final Collection<IBDState> overlappingIBDRegions = ibdStateMaps.get(siblingPair).getOverlapping(interval);

                    processSingleInterval(vc, siblingPair, gi, interval, overlappingIBDRegions, paralogId, minSibLL, minFamilyLL);

                } else {
                    throw new GATKException("No m2 mode specified");
                }
            }
        }
        return null;
    }

    public void close() {
        if (endsConfidenceWriter != null) {
            endsConfidenceWriter.close();
        }
        if (reportWriter != null) {
            reportWriter.close();
        }
    }

    private LikelihoodCombination getBestLikelihood(final VariantContext vc, final VariantContext end1Vc, final VariantContext end2Vc, final String sample) {
        final Genotype sib1Gt = vc.getGenotype(sample);
        final double[] sib1SDLikelihoods = GenotypeUtilities.getCopyNumberLogLikelihoods(sib1Gt);
        final double[] sib1End1Likelihoods = GenotypeUtilities.getCopyNumberLogLikelihoods(end1Vc.getGenotype(sample));
        final double[] sib1End2Likelihoods = GenotypeUtilities.getCopyNumberLogLikelihoods(end2Vc.getGenotype(sample));
        if (sib1End1Likelihoods == null || sib1End2Likelihoods == null) {
            return null;
        }

        final List<LikelihoodCombination> likelihoodCombinations = new ArrayList<>(sib1End1Likelihoods.length * sib1End2Likelihoods.length * sib1SDLikelihoods.length);
        //double totalConsistentLL = Math.log10(0);
        for (int i = 0; i < sib1End1Likelihoods.length; i++) {
            for (int j = 0; j < sib1End2Likelihoods.length; j++) {
                for (int k = 0; k < sib1SDLikelihoods.length; k++) {
                    final double logLikelihood = sib1End1Likelihoods[i] + sib1End2Likelihoods[j] + sib1SDLikelihoods[k];
                    if (i + j == k) {
                        likelihoodCombinations.add(new LikelihoodCombination(i, j, k, logLikelihood));
                        //totalConsistentLL = MathUtils.log10sumLog10(new double[] {totalConsistentLL, logLikelihood});
                    }
                }
            }
        }
        Collections.sort(likelihoodCombinations, new Comparator<LikelihoodCombination>() {
            @Override
            public int compare(final LikelihoodCombination o1, final LikelihoodCombination o2) {
                return new Double(o1.ll).compareTo(o2.ll);
            }
        });

        if (likelihoodCombinations.size() < 2) {
            return null;
        }
        final LikelihoodCombination bestLLCombo = likelihoodCombinations.get(likelihoodCombinations.size() - 1);
        final LikelihoodCombination nextBestLLCombo = likelihoodCombinations.get(likelihoodCombinations.size() - 2);
        final LikelihoodCombination bestLLComboNormalized = new LikelihoodCombination(bestLLCombo.end1Cn,
                bestLLCombo.end2Cn,
                bestLLCombo.sdCn,
                bestLLCombo.ll - nextBestLLCombo.ll);
        return bestLLComboNormalized;
    }

    private void processDoubleIntervals(final VariantContext vc, final SiblingPair siblingPair, final GenomeInterval gi) {
        final Interval interval;
        String dupIntervalString = vc.getAttributeAsString("DUPINTERVALS", null);


        interval = new Interval(gi.getSequenceName(), gi.getStart(), gi.getEnd());
        GenomeInterval dupGi = GenomeInterval.parse(dupIntervalString);
        Interval dupInterval = new Interval(dupGi.getSequenceName(), dupGi.getStart(), dupGi.getEnd());


        final Collection<IBDState> overlappingIBDRegions = ibdStateMaps.get(siblingPair).getOverlapping(interval);
        final Collection<IBDState> dupOverlappingIBDRegions = ibdStateMaps.get(siblingPair).getOverlapping(dupInterval);

        if (overlappingIBDRegions.size() == 1 && dupOverlappingIBDRegions.size() == 1) {
            final Interval overlappingIBDInterval = ibdStateMaps.get(siblingPair).getOverlappingKeys(interval).iterator().next();
            final IBDState ibdState = overlappingIBDRegions.iterator().next();

            final Interval dupOverlappingIBDInterval = ibdStateMaps.get(siblingPair).getOverlappingKeys(dupInterval).iterator().next();
            final IBDState dupIbdState = dupOverlappingIBDRegions.iterator().next();

            if (overlappingIBDInterval.getIntersectionLength(interval) < interval.length()) {
                System.out.println("IBD region " + overlappingIBDInterval + ", " + ibdState + " does not fully contain " + vc.getID() + " in sib pair " + siblingPair);
                writeReportLine(vc, siblingPair, null, null, "PARTIAL_OVERLAP");
                return;
            }

            if (dupOverlappingIBDInterval.getIntersectionLength(dupInterval) < dupInterval.length()) {
                System.out.println("Dup IBD region " + dupOverlappingIBDInterval + ", " + dupIbdState + " does not fully contain DUPINTERVAL " + dupInterval + " of " + vc.getID() + " in sib pair " + siblingPair);
                writeReportLine(vc, siblingPair, null, null, "PARTIAL_OVERLAP_DUPINTERVAL");
                return;
            }

            if (!ibdState.getRefinedState().equals(dupIbdState.getRefinedState())) {
                writeReportLine(vc, siblingPair, null, null, "DIFF_DUPINTERVAL_IBD_STATES");
                return;
            }

            final Integer distanceToIBDStateChange = Math.min(
                    Math.min(interval.getStart() - overlappingIBDInterval.getStart(), overlappingIBDInterval.getEnd() - interval.getEnd()),
                    Math.min(dupInterval.getStart() - dupOverlappingIBDInterval.getStart(), dupOverlappingIBDInterval.getEnd() - dupInterval.getEnd())
            );
            writeReportLine(vc, siblingPair, ibdState, distanceToIBDStateChange, null);
        } else {
            if (overlappingIBDRegions.size() > 1) {
                System.out.println("Multiple overlapping IBD regions for variant " + vc.getID() + " in sib pair " + siblingPair + ": " + ibdStateMaps.get(siblingPair).getOverlappingKeys(interval));
                writeReportLine(vc, siblingPair, null, null, "MULTIPLE_IBD_REGIONS");
                return;
            }
            if (dupOverlappingIBDRegions.size() > 1) {
                System.out.println("Multiple overlapping IBD regions for DUPINTERVAL of variant " + vc.getID() + " in sib pair " + siblingPair + ": " + ibdStateMaps.get(siblingPair).getOverlappingKeys(dupInterval));
                writeReportLine(vc, siblingPair, null, null, "MULTIPLE_DUPINTERVAL_IBD_REGIONS");
                return;
            }
            if (overlappingIBDRegions.size() == 0) {
                writeReportLine(vc, siblingPair, null, null, "NO_IBD_REGIONS");
            }
            if (dupOverlappingIBDRegions.size() == 0) {
                writeReportLine(vc, siblingPair, null, null, "NO_DUPINTERVAL_IBD_REGIONS");
            }

        }
    }

    private void processSingleInterval(final VariantContext vc, final SiblingPair siblingPair, final GenomeInterval gi,
                                       final Interval interval, final Collection<IBDState> overlappingIBDRegions, final String intervalId, final Double minSibLL, final Double minFamilyLL) {
        if (overlappingIBDRegions.size() == 1) {
            final Interval overlappingIBDInterval = ibdStateMaps.get(siblingPair).getOverlappingKeys(interval).iterator().next();
            final IBDState ibdState = overlappingIBDRegions.iterator().next();
            String errorCond = null;
            if (minSibLL != null) {
                if (ibdState.getRefinedState().equals("TWO")) {
                    if (minSibLL < endConfidenceLikelihoodThreshold) {
                        errorCond = "IBD2_AND_MIN_SIB_PARA_LIKELIHOOD_THRESHOLD";
                    }
                } else {
                    if (minFamilyLL < endConfidenceLikelihoodThreshold) {
                        errorCond = "MIN_PARENT_PARA_LIKELIHOOD_THRESHOLD";
                    }
                }
            }

            if (overlappingIBDInterval.getIntersectionLength(interval) < interval.getEnd() - interval.getStart()) {
                final Interval preSiteInterval = new Interval(gi.getSequenceName(), 1, gi.getStart() - 1);
                if (ibdStateMaps.get(siblingPair).getOverlapping(preSiteInterval).isEmpty()) {
                    writeReportLine(vc, siblingPair, ibdState, null, errorCond);
                } else {
                    System.out.println("IBD region " + overlappingIBDInterval + ", " + ibdState + " does not fully contain " + intervalId + " in sib pair " + siblingPair);
                    writeReportLine(vc, siblingPair, null, null, "PARTIAL_OVERLAP");
                    return;
                }
            }
            final Integer distanceToIBDStateChange = Math.min(interval.getStart() - overlappingIBDInterval.getStart(), overlappingIBDInterval.getEnd() - interval.getEnd());

            writeReportLine(vc, siblingPair, ibdState, distanceToIBDStateChange, errorCond);
        } else {
            if (overlappingIBDRegions.size() > 1) {
                System.out.println("Multiple overlapping IBD regions for interval " + intervalId + " in sib pair " + siblingPair + ": " + ibdStateMaps.get(siblingPair).getOverlappingKeys(interval));
                writeReportLine(vc, siblingPair, null, null, "MULTIPLE_IBD_REGIONS");
                return;
            }
            // check to see if this is an uncovered region at the telomere, in which case try targeting it using the next
            // IBD region on the chromosome
            final Interval preSiteInterval = new Interval(gi.getSequenceName(), 1, gi.getStart() - 1);
            if (ibdStateMaps.get(siblingPair).getOverlapping(preSiteInterval).isEmpty()) {
                // use getDebugTree to get the raw interval tree and its min
                final IntervalTree<IBDState> treeForContig = ibdStateMaps.get(siblingPair).debugGetTree(gi.getSequenceName());
                if (treeForContig != null) {
                    final IntervalTree.Node<IBDState> minNode = treeForContig.min();
                    final IBDState nextIBDState = minNode.getValue();
                    String errorCond = null;
                    if (minSibLL != null) {
                        if (nextIBDState.getRefinedState().equals("TWO")) {
                            if (minSibLL < endConfidenceLikelihoodThreshold) {
                                errorCond = "IBD2_AND_MIN_SIB_PARA_LIKELIHOOD_LT_THRESHOLD";
                            }
                        } else {
                            if (minFamilyLL < endConfidenceLikelihoodThreshold) {
                                errorCond = "MIN_PARENT_PARA_LIKELIHOOD_LT_THRESHOLD";
                            }
                        }
                    }
                    writeReportLine(vc, siblingPair, nextIBDState, null, errorCond);
                    return;
                }
            }

            writeReportLine(vc, siblingPair, null, null, "NO_IBD_REGIONS");
        }
    }

    private boolean isM2Genotyping(final VariantContext vc) {
        return vc.getID().startsWith("GS_SD_M2");
    }

    public static class SiblingPair {
        final String sib1;
        final String sib2;

        public SiblingPair(final String sib1, final String sib2) {
            if (sib1.compareTo(sib2) < 1) {
                this.sib1 = sib1;
                this.sib2 = sib2;
            } else {
                this.sib2 = sib1;
                this.sib1 = sib2;
            }
        }

        public String getName() {
            return sib1 + "-" + sib2;
        }

        @Override
        public String toString() {
            return getName();
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final SiblingPair that = (SiblingPair) o;

            if (sib1 != null ? !sib1.equals(that.sib1) : that.sib1 != null) return false;
            return !(sib2 != null ? !sib2.equals(that.sib2) : that.sib2 != null);

        }

        @Override
        public int hashCode() {
            int result = sib1 != null ? sib1.hashCode() : 0;
            result = 31 * result + (sib2 != null ? sib2.hashCode() : 0);
            return result;
        }
    }

    // todo: somehow unify this code with what's in the SiblingIBD annotation package in the GATK codebase
    private static class IBDState {
        private String state;
        private Integer sites;
        private Double paternalIndicators;
        private Double maternalIndicators;

        public IBDState(final String state, final Integer sites, Double paternalIndicators, final Double maternalIndicators) {
            this.state = state;
            this.sites = sites;
            this.paternalIndicators = paternalIndicators;
            this.maternalIndicators = maternalIndicators;
        }

        @Override
        public boolean equals(final Object o) {

            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final IBDState ibdState = (IBDState) o;

            if (state != null ? !state.equals(ibdState.state) : ibdState.state != null) return false;
            if (sites != null ? !sites.equals(ibdState.sites) : ibdState.sites != null)
                return false;

            if (paternalIndicators != null ? !paternalIndicators.equals(ibdState.paternalIndicators) : ibdState.paternalIndicators != null)
                return false;
            return !(maternalIndicators != null ? !maternalIndicators.equals(ibdState.maternalIndicators) : ibdState.maternalIndicators != null);

        }

        @Override
        public int hashCode() {
            int result = state != null ? state.hashCode() : 0;
            result = 31 * result + (sites != null ? sites.hashCode() : 0);
            result = 31 * result + (paternalIndicators != null ? paternalIndicators.hashCode() : 0);
            result = 31 * result + (maternalIndicators != null ? maternalIndicators.hashCode() : 0);
            return result;
        }

        public String getRefinedState() {
            if (! "ONE".equals(state)) {
                return state;
            } else {
                if (paternalIndicators > maternalIndicators) {
                    return "ONE-P";
                } else if (maternalIndicators > paternalIndicators) {
                    return "ONE-M";
                } else {
                    return "ONE-?";
                }
            }

        }

        public String getParentConfidence() {
            if (! "ONE".equals(state)) {
                return "NA";
            } else {
                if (paternalIndicators > maternalIndicators) {
                    return String.valueOf(paternalIndicators / (paternalIndicators + maternalIndicators));
                } else if (maternalIndicators > paternalIndicators) {
                    return String.valueOf(maternalIndicators / (paternalIndicators + maternalIndicators));
                } else {
                    return "NA";
                }
            }

        }

        public Integer getSites() {
            return sites;
        }

        @Override
        public String toString() {
            return "IBDState{" +
                    "state='" + state + '\'' +
                    ", sites=" + sites +
                    ", paternalIndicators=" + paternalIndicators +
                    ", maternalIndicators=" + maternalIndicators +
                    '}';
        }
    }

    private class LikelihoodCombination {

        public LikelihoodCombination(final int end1Cn, final int end2Cn, final int sdCn, final double ll) {
            this.end1Cn = end1Cn;
            this.end2Cn = end2Cn;
            this.sdCn = sdCn;
            this.ll = ll;
        }

        public int getEnd1Cn() {
            return end1Cn;
        }

        public int getEnd2Cn() {
            return end2Cn;
        }

        public int getSdCn() {
            return sdCn;
        }

        public double getLl() {
            return ll;
        }

        final int end1Cn;
        final int end2Cn;
        final int sdCn;
        final double ll;

    }
}