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


import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.samples.MendelianViolation;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by cwhelan on 5/27/14.
 *
 * Tracks observations of quartet genotypes across the genome and executes the HMM model when
 * appropriate.
 *
 */
public class QuartetIBDStateModel implements IBDStateModel {

    private final Map<SiblingPair, QuartetIBDStateHMM> quartetIBDStateTrackers = new HashMap<>();
    private final PrintStream spFile;
    private final Integer gqThreshold;
    private final Double hetABThreshold;
    private String chrom = "NA";

    public QuartetIBDStateModel(final List<SiblingPair> siblingPairs,
                                final Integer gqThreshold,
                                final Double hetABThreshold,
                                final PrintStream spFile) {
        for (final SiblingPair siblingPair : siblingPairs) {
            quartetIBDStateTrackers.put(siblingPair, new QuartetIBDStateHMM());
        }
        this.gqThreshold = gqThreshold;
        this.hetABThreshold = hetABThreshold;
        this.spFile = spFile;
    }

    @Override
    public List<VariantIBDState> addSite(final VariantContext vc, final GenotypesContext genotypes, final SiblingPair siblingPair, final Genotype sib1Gt, final Genotype sib2Gt, final List<VariantContext> overlappingCNVs, final boolean keepVC) {
        final List<VariantIBDState> result = new ArrayList<>();

        if (! "NA".equals(chrom) && ! vc.getContig().equals(chrom)) {
            for (SiblingPair modeledPair : quartetIBDStateTrackers.keySet()) {
                final QuartetIBDStateHMM ibdStateHmm = quartetIBDStateTrackers.get(modeledPair);
                final Iterator<IBDLocus> iterator = ibdStateHmm.runModelAndCoallateResults();

                while (iterator.hasNext()) {
                    final IBDLocus locus = iterator.next();
                    addResult(modeledPair, result, locus, keepVC);
                }
            }
        }
        chrom = vc.getContig();

        final String father = siblingPair.sib1.getPaternalID();
        final String mother = siblingPair.sib1.getMaternalID();

        if (genotypes.containsSample(father) && genotypes.containsSample(mother)) {
            final Genotype paternalGt = genotypes.get(father);
            final Genotype maternalGt = genotypes.get(mother);
            if (MendelianViolation.isViolation(paternalGt, maternalGt, sib1Gt) ||
                    MendelianViolation.isViolation(paternalGt, maternalGt, sib2Gt)) {
                return result;
            }

            if (passesGQThreshold(sib1Gt, sib2Gt, paternalGt, maternalGt) &&
                passesABHetThreshold(sib1Gt, sib2Gt, paternalGt, maternalGt) &&
                noOverlappingCNV(overlappingCNVs,
                        siblingPair.sib1.getID(),
                        siblingPair.sib2.getID(),
                        father,
                        mother)) {
                final QuartetIBDStateHMM ibdStateHmm = quartetIBDStateTrackers.get(siblingPair);
                final IBDObservation observation = IBDObservation.getIBDStateObservation(paternalGt, maternalGt, sib1Gt, sib2Gt);
                if (observation == IBDObservation.ZERO_OR_ONE_OR_TWO) {
                    return result;
                }

                if (spFile != null) {
                    spFile.print(vc.getContig() + "\t" + vc.getStart() + "\t" + siblingPair.getName() + "\t" + observation.ordinal() + "\n");
                }

                ibdStateHmm.addObservation(vc, paternalGt, maternalGt, observation, getParentalIndicator(siblingPair, vc), keepVC);

            }
        }
        return result;
    }

    private boolean noOverlappingCNV(final List<VariantContext> overlappingCNVs, final String ... samples) {
        for (final VariantContext cnvVc : overlappingCNVs) {
            if (cnvVc.isFiltered()) {
                continue;
            }
            for (final String sample : samples) {
                if (cnvVc.hasGenotype(sample)) {
                    final Genotype gt = cnvVc.getGenotype(sample);
                    if (gt.hasAnyAttribute("CN")) {
                        final Integer cn = Integer.valueOf((String) gt.getAnyAttribute("CN"));
                        if (cn != 2) {
                            return false;
                        }
                    }
                }
            }
        }
        return true;
    }

    private boolean passesABHetThreshold(final Genotype sib1Gt, final Genotype sib2Gt, final Genotype paternalGt, final Genotype maternalGt) {
        return (! sib1Gt.isHet() || passesABHetThreshold(sib1Gt)) &&
                (! sib2Gt.isHet() || passesABHetThreshold(sib2Gt)) &&
                (! paternalGt.isHet() || passesABHetThreshold(paternalGt)) &&
                (! maternalGt.isHet() || passesABHetThreshold(maternalGt));
    }

    private boolean passesABHetThreshold(final Genotype gt) {
        final double rawAB = ((double) gt.getAD()[0]) / (gt.getAD()[0] + gt.getAD()[1]);
        final double AB = Math.min(rawAB, 1 - rawAB);
        return AB > hetABThreshold;
    }

    private boolean passesGQThreshold(final Genotype sib1Gt, final Genotype sib2Gt, final Genotype paternalGt, final Genotype maternalGt) {
        return gqThreshold <= Math.min(Math.min(sib1Gt.getGQ(), sib2Gt.getGQ()), Math.min(paternalGt.getGQ(), maternalGt.getGQ()));
    }

    private void addResult(final SiblingPair siblingPair, final List<VariantIBDState> result, final IBDLocus locus, final boolean keepVc) {

        result.add(new VariantIBDState(locus.chr, locus.pos,
                keepVc ? locus.vc : null,
                siblingPair, locus.state, locus.statePosteriors, locus.parentalIndicator,
                locus.observationAgrees, locus.observation));
    }

    private IBDParentalAgreement getParentalIndicator(final SiblingPair siblingPair, final VariantContext vc) {
        final String father = siblingPair.sib1.getPaternalID();
        final String mother = siblingPair.sib1.getMaternalID();

        final Genotype paternalGt = vc.getGenotype(father);
        final Genotype maternalGt = vc.getGenotype(mother);

        final Genotype sib1Gt = vc.getGenotype(siblingPair.sib1.getID());
        final Genotype sib2Gt = vc.getGenotype(siblingPair.sib2.getID());


        return getParentalIndicator(paternalGt, maternalGt, sib1Gt, sib2Gt);
    }

    /**
     * Assumes IBD1. Returns -1 if parental shared IBD can not be determined, 0 if sibs share paternal allele, 1 if sibs share
     * maternal allele.
     * @param paternalGt
     * @param maternalGt
     * @param sib1Gt
     * @param sib2Gt
     * @return
     */
    private IBDParentalAgreement getParentalIndicator(final Genotype paternalGt, final Genotype maternalGt, final Genotype sib1Gt, final Genotype sib2Gt) {
        if (maternalGt.isHet() && paternalGt.isHom()) {
            if (sib1Gt.isHet() && sib2Gt.isHet()) {
                return IBDParentalAgreement.MATERNAL;
            }
            if (sib1Gt.isHet() && !sib2Gt.isHet() || !sib1Gt.isHet() && sib2Gt.isHet()) {
                return IBDParentalAgreement.PATERNAL;
            }
        }
        if (paternalGt.isHet() && maternalGt.isHom()) {
            if (sib1Gt.isHet() && sib2Gt.isHet()) {
                return IBDParentalAgreement.PATERNAL;
            }
            if (sib1Gt.isHet() && !sib2Gt.isHet() || !sib1Gt.isHet() && sib2Gt.isHet()) {
                return IBDParentalAgreement.MATERNAL;
            }
        }

        return IBDParentalAgreement.UNKNOWN;
    }

    @Override
    public List<VariantIBDState> finalizeModel(final boolean keepVc) {
        final List<VariantIBDState> finalValues = new ArrayList<>();
        for (final SiblingPair siblingPair : quartetIBDStateTrackers.keySet()) {
            final QuartetIBDStateHMM quartetIBDStateHMM = quartetIBDStateTrackers.get(siblingPair);
            final Iterator<IBDLocus> iterator = quartetIBDStateHMM.runLastChrom();
            while (iterator.hasNext()) {
                final IBDLocus locus = iterator.next();
                addResult(siblingPair, finalValues, locus, keepVc);
            }
        }
        return finalValues;
    }
}
