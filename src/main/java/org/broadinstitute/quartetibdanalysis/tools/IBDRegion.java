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
 * Represents a genome interval and its IBD state
 */
public class IBDRegion {
    String chr;
    int start;
    int end;
    IBDState state;
    int paternalIndicators;
    int maternalIndicators;
    int sites;
    int siteAgreements;

    public IBDRegion(final String chr, final int start, final int end, final IBDState state, final int paternalIndicators, final int maternalIndicators, final int sites, final int siteAgreements) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.state = state;
        this.paternalIndicators = paternalIndicators;
        this.maternalIndicators = maternalIndicators;
        this.sites = sites;
        this.siteAgreements = siteAgreements;
    }


    @Override
    public String toString() {
        return "IBDRegion{" +
                "chr='" + chr + '\'' +
                ", start=" + start +
                ", end=" + end +
                ", state=" + state +
                ", paternalIndicators=" + paternalIndicators +
                ", maternalIndicators=" + maternalIndicators +
                ", sites=" + sites +
                ", siteAgreements=" + siteAgreements +
                '}';
    }
}
