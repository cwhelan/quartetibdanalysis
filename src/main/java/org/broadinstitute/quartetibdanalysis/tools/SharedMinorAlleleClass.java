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
 * Class to represent the state of unambiguous minor allele sharing at a given locus between two siblings.
 *
*/
public enum SharedMinorAlleleClass {
    HOMVAR_HOMVAR,
    HETVAR_HETVAR,
    HOMVAR_HETVAR,
    HETVAR_HOMREF,
    HOMREF_HOMREF,
    HOMVAR_HOMREF
}
