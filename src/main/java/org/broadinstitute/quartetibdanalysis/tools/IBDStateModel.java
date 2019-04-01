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

import java.util.List;

/**
 * Created by cwhelan on 5/27/14.
 *
 * An interface to support varying models for predicting IBD state across the genome.
 *
 */
public interface IBDStateModel {

    /**
     * Add a genotyped site to the model.
     *  @param vc
     * @param genotypes
     * @param siblingPair
     * @param sib1Gt
     * @param sib2Gt    @return an empty list if we can't make any predictions right now; otherwise a list of sites at which we are ready to make IBD predictions
     * @param overlappingCNVs
     * */
    public List<VariantIBDState> addSite(final VariantContext vc, final GenotypesContext genotypes, final SiblingPair siblingPair, final Genotype sib1Gt, final Genotype sib2Gt, final List<VariantContext> overlappingCNVs, final boolean keepVC);

    /**
     * Call after all sites have been processed
     *
     * @return a list of any additional sites for which we can make IBD predictions
     */
    public List<VariantIBDState> finalizeModel(final boolean keepVc);
}
