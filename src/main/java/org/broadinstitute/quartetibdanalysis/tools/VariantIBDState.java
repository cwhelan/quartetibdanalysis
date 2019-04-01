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
* Created by cwhelan on 5/27/14.
*
 * The IBD state at locus for a particular sibling pair
 *
*/
public class VariantIBDState {
    final String chr;
    final int pos;
    final VariantContext vc;
    final SiblingPair siblingPair;
    final IBDState ibdState;
    final double[] statePosteriors;
    final IBDParentalAgreement parentalAgreement;
    final boolean agrees;
    final IBDObservation observation;

    public VariantIBDState(final String chr, final int pos, final VariantContext vc, final SiblingPair siblingPair,
                           final IBDState ibdState, final double[] statePosteriors, final IBDParentalAgreement parentalAgreement,
                           final boolean agrees, final IBDObservation observation) {
        this.chr = chr;
        this.pos = pos;
        this.vc = vc;
        this.statePosteriors = statePosteriors;
        this.ibdState = ibdState;
        this.siblingPair = siblingPair;
        this.parentalAgreement = parentalAgreement;
        this.agrees = agrees;
        this.observation = observation;
    }

    @Override
    public String toString() {
        return "VariantIBDState: " + chr + "\t" + pos + "\t" + siblingPair + "\t" + ibdState;
    }
}
