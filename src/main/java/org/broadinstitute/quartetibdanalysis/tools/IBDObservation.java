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
import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * Created by cwhelan on 5/8/14.
 *
 * Given the genotypes of the two parents and two siblings at a given locus,
 * this class represents the type of information about IBD state that can be
 * inferred. Based on the genotypes, each site can imply IBD0, IBD1, or IBD2,
 * or one of two of the IBD states, or leave the possibility of any of the three
 * IBD states.
 */
public enum IBDObservation {
    ZERO,
    ONE,
    TWO,
    ZERO_OR_ONEF,
    ZERO_OR_ONEM,
    ZERO_OR_TWO,
    ONEF_OR_TWO,
    ONEM_OR_TWO,
    ZERO_OR_ONE_OR_TWO;

    private static boolean bothHom(final Genotype sib1Gt, final Genotype sib2Gt) {
        return sib1Gt.isHom() && sib2Gt.isHom();
    }

    private static boolean bothHet(final Genotype sib1Gt, final Genotype sib2Gt) {
        return sib1Gt.isHet() && sib2Gt.isHet();
    }

    /**
     * Given the four genotypes at a locus, return the type of IBD state observation
     *
     * @param paternalGt paternal genotype
     * @param maternalGt maternal genotype
     * @param s1 sibling genotype 1
     * @param s2 sibling genotype 2
     * @return the type of IBD observation indicated by the genotype configuration of the quartet
     *
     */
    public static IBDObservation getIBDStateObservation(final Genotype paternalGt, final Genotype maternalGt, final Genotype s1, final Genotype s2) {
        if (bothHet(paternalGt, maternalGt)) {
            if (bothHom(s1, s2)) {
                if (s1.sameGenotype(s2)) {
                    return IBDObservation.TWO;
                } else {
                    return IBDObservation.ZERO;
                }
            } else if (bothHet(s1,s2)) {
                return IBDObservation.ZERO_OR_TWO;
            } else {
                // one sib het, the other hom
                return IBDObservation.ONE;
            }
        } else if (paternalGt.isHet() && maternalGt.isHom() || (paternalGt.isHom() && maternalGt.isHet())) {
            if (bothHet(s1, s2) || bothHom(s1, s2)) {
                if (paternalGt.isHet()) {
                    return IBDObservation.ONEF_OR_TWO;
                } else {
                    return IBDObservation.ONEM_OR_TWO;
                }
            } else {
                // one sib het, the other hom
                if (paternalGt.isHom()) {
                    return IBDObservation.ZERO_OR_ONEF;
                } else {
                    return IBDObservation.ZERO_OR_ONEM;
                }
            }
        } else {
            if (! bothHom(paternalGt,maternalGt)) {
                throw new GATKException("Unexpected parental genotypes " + paternalGt + " and " + maternalGt + "; expect both to be het, both to be hom, or one het and one hom");
            }
            assert(bothHom(paternalGt,maternalGt));
            return IBDObservation.ZERO_OR_ONE_OR_TWO;
        }
    }

    public boolean agreesWith(IBDState state) {
        switch(this) {
            case ZERO: return state == IBDState.ZERO;
            case ONE: return state == IBDState.ONEF || state == IBDState.ONEM;
            case TWO: return state == IBDState.TWO;
            case ZERO_OR_ONEF:  return state == IBDState.ZERO || state == IBDState.ONEF;
            case ZERO_OR_ONEM:  return state == IBDState.ZERO || state == IBDState.ONEM;
            case ONEF_OR_TWO:  return state == IBDState.ONEF || state == IBDState.TWO;
            case ONEM_OR_TWO:  return state == IBDState.ONEM || state == IBDState.TWO;
            case ZERO_OR_TWO: return state == IBDState.ZERO || state == IBDState.TWO;
            case ZERO_OR_ONE_OR_TWO: return true;
        }
        return false;
    }
}
