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


import org.broadinstitute.hellbender.utils.samples.Sample;

/**
* Created by cwhelan on 5/27/14.
 *
 * Organizes a pair of siblings for IBD analysis.
 *
*/
class SiblingPair {
    final Sample sib1;
    final Sample sib2;

    public SiblingPair(final Sample sib1, final Sample sib2) {
        this.sib1 = sib1;
        this.sib2 = sib2;
    }

    public String getName() {
        return sib1.getID() + "-" + sib2.getID();
    }

    @Override
    public String toString() {
        return getName();
    }
}
