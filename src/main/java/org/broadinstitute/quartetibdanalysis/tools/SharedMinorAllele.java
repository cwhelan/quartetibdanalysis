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
* Created by cwhelan on 5/1/14.
 *
 * Represents the number of minor alleles unambiguously shared between two siblings at a locus,
 * as well as the total number of minor alleles at that locus.
 *
*/
public class SharedMinorAllele {
    int matches;
    int max;

    public SharedMinorAllele(final int matches, final int max) {
        this.matches = matches;
        this.max = max;
    }

    public SharedMinorAlleleClass sharedAlleleClass() {
        if (max == 2) {
            if (matches == 2) {
                return SharedMinorAlleleClass.HOMVAR_HOMVAR;
            } else if (matches == 1) {
                return SharedMinorAlleleClass.HOMVAR_HETVAR;
            } else {
                return SharedMinorAlleleClass.HOMVAR_HOMREF;
            }
        } else if (max == 1) {
            if (matches == 1) {
                return SharedMinorAlleleClass.HETVAR_HETVAR;
            } else {
                return SharedMinorAlleleClass.HETVAR_HOMREF;
            }
        } else {
            return SharedMinorAlleleClass.HOMREF_HOMREF;
        }
    }
}
