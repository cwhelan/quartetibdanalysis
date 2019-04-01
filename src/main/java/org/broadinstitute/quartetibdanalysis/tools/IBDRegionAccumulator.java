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

/**
 * Created by cwhelan on 5/6/14.
 *
 * When fed observations loci with their predicted IBD states, emits larger regions
 * across which all IBD predictions are consistent.
 *
 */
public class IBDRegionAccumulator {

    String chr;
    int start;
    int currentMax;
    IBDState state;
    int paternalCount;
    int maternalCount;
    int sites;
    int siteAgreements;

    /**
     * Given a locus and IBD state, have we changed state from the previous set of observations?
     *
     * @param chr
     * @param coord
     * @param state
     * @param agrees
     * @return null if we are still in the same state; an IBDRegion if we have changed state
     */
    public IBDRegion regionChange(final String chr, final int coord, final IBDState state, final IBDParentalAgreement parentalAgreement, final boolean agrees) {
        if (this.chr == null) {
            resetState(chr, coord, state);
            if (parentalAgreement == IBDParentalAgreement.PATERNAL) {
                this.paternalCount++;
            }
            if (parentalAgreement == IBDParentalAgreement.MATERNAL) {
                this.maternalCount++;
            }
            if (agrees) {
                this.siteAgreements++;
            }
            this.sites++;
            return null;
        }
        final boolean regionChanged = (!chr.equals(this.chr) || state != this.state);
        IBDRegion result = null;
        if (regionChanged) {
            result = new IBDRegion(this.chr, this.start, this.currentMax, this.state, this.paternalCount, this.maternalCount, this.sites, this.siteAgreements);
            resetState(chr, coord, state);
        }
        this.chr = chr;
        this.currentMax = coord;
        this.sites++;
        if (parentalAgreement == IBDParentalAgreement.PATERNAL) {
            this.paternalCount++;
        }
        if (parentalAgreement == IBDParentalAgreement.MATERNAL) {
            this.maternalCount++;
        }
        if (agrees) {
            this.siteAgreements++;
        }

        return result;
    }

    /**
     * @return the final IBD region in the genome, when there are no more loci to accumulate
     */
    public IBDRegion getFinalRegion() {
        return new IBDRegion(this.chr, this.start, this.currentMax, this.state, this.paternalCount, this.maternalCount, this.sites, this.siteAgreements);
    }

    private void resetState(final String chr, final int coord, final IBDState state) {
        this.chr = chr;
        this.start = coord;
        this.currentMax = coord;
        this.state = state;
        this.paternalCount = 0;
        this.maternalCount = 0;
        this.sites = 0;
        this.siteAgreements = 0;
    }

}
