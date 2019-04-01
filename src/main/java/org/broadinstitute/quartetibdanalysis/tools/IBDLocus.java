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

import htsjdk.variant.variantcontext.VariantContext;

/**
* Created by cwhelan on 5/20/14.
 *
 * A simple data structure to hold the the predicted IBD state at a given locus.
 *
*/
public class IBDLocus {
    String chr;
    int pos;
    VariantContext vc;
    IBDState state;
    double[] statePosteriors;
    IBDParentalAgreement parentalIndicator;
    boolean observationAgrees;
    IBDObservation observation;

    public IBDLocus(final String chr, final int pos, final VariantContext vc,
                    final IBDState state, final double[] statePosteriors, final IBDParentalAgreement parentalIndicator,
                    final boolean observationAgrees, final IBDObservation observation) {
        this.chr = chr;
        this.pos = pos;
        this.vc = vc;
        this.state = state;
        this.statePosteriors = statePosteriors;
        this.parentalIndicator = parentalIndicator;
        this.observationAgrees = observationAgrees;
        this.observation = observation;
    }
}
