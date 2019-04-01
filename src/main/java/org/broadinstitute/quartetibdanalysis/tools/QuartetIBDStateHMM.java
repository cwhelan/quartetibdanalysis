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
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.*;

/**
 * Created by cwhelan on 5/19/14.
 *
 * Models IBD state across the genomes of two siblings using an HMM.
 *
 */
public class QuartetIBDStateHMM {
    final protected static Logger logger = Logger.getLogger(QuartetIBDStateHMM.class);

    private final List<IBDObservation> observations = new ArrayList<>();
    private final List<VariantContext> sites = new ArrayList<>();
    private final List<String> seqs = new ArrayList<>();
    private final List<Integer> starts = new ArrayList<>();
    private final List<IBDParentalAgreement> parentalIndicators = new ArrayList<>();

    // these constants tuned to be arbitrarily low based on experiments conducted on
    // a human whole genome data set; we found that the system performs best when
    // it is discouraged from switching from one state to another.
    private static final double DEFAULT_STATE_CHANGE_PROB = 1.0 / 10000000000000000.0;
    private static final double DEFAULT_ERROR_PROB = 1.0 / 1000000000000000000.0;

    // emission probabilities for the HMM; we found that uniform emission probabilities
    // across valid or invalid observations performed well on a human test data set,
    // with a slight tweak to give extra probability to IBD_ZERO_OR_TWO observations in
    // IBD1 states - a common error mode due to overcalling of het genotypes
    public static final double NON_ERROR_EMISSION_PROBABILITY = .97;
    public static final double ERROR_EMISSION_PROBABILITY = .03;
    public static final double IBD1_HET_ERROR_EMISSION_PROBABILITY = .02;
    public static final double IDB1_OTHER_ERROR_EMISSION_PROBABILITY = .01;

    private String chrom = null;

    private EnumMap<IBDState, Double> initialProbs;
    private EnumMap<IBDState, Double[]> transitionProbs;

    public QuartetIBDStateHMM() {
        initProbs(DEFAULT_STATE_CHANGE_PROB, DEFAULT_ERROR_PROB);
    }

    public QuartetIBDStateHMM(final double recombinationProb, final double errorProb) {
        initProbs(recombinationProb, errorProb);
    }

    /**
     * Initialize initial and transition probabilities
     *
     * HMM probabilities are modeled simplistically, with uniform initial probabilities and
     * uniform probabilities for switching states
     * @param recombinationProb
     * @param errorProb
     */
    private void initProbs(final double recombinationProb, final double errorProb) {

        initialProbs = new EnumMap<>(IBDState.class);
        initialProbs.put(IBDState.ZERO, Math.log10(.33));
        initialProbs.put(IBDState.ONE, Math.log10(.33));
        initialProbs.put(IBDState.TWO, Math.log10(.33));

        transitionProbs = new EnumMap<>(IBDState.class);
        transitionProbs.put(IBDState.ZERO, new Double[]{Math.log10((1 - recombinationProb) - (errorProb / 2)), Math.log10(recombinationProb - (errorProb / 2)), Math.log10(errorProb)});
        transitionProbs.put(IBDState.ONE, new Double[]{Math.log10((recombinationProb / 2)), Math.log10(1 - recombinationProb), Math.log10(recombinationProb / 2)});
        transitionProbs.put(IBDState.TWO, new Double[]{Math.log10(errorProb), Math.log10(recombinationProb - (errorProb / 2)), Math.log10((1 - recombinationProb) - (errorProb / 2))});
    }

    /**
     * Add an observation to the HMM
     *
     * @param vc
     * @param observation
     * @return
     */
    public Iterator<IBDLocus> addObservation(final VariantContext vc, final IBDObservation observation, final IBDParentalAgreement parentalIndicator, boolean keepVc) {
        observations.add(observation);

        if (keepVc) {
            sites.add(vc);
        }

        parentalIndicators.add(parentalIndicator);
        seqs.add(vc.getContig());
        starts.add(vc.getStart());

        return Collections.<IBDLocus>emptyList().iterator();
    }

    /**
     * Execute Viterbi for the final chromosome and return the results
     * @return
     */
    public Iterator<IBDLocus> runLastChrom() {
        return runModelAndCoallateResults();
    }

    /**
     * execute Viterbi
     * and return an iterator over the predicted IBD state for each locus for which we've seen an observation. Otherwise
     */
    public Iterator<IBDLocus> runModelAndCoallateResults() {
        final ModelResults modelResults = runModel();
        final IBDState[] states = modelResults.viterbiStates;
        final List<IBDLocus> results = new ArrayList<>(observations.size());
        for (int i = 0; i < states.length; i++) {
            VariantContext vc = null;
            if (sites.size() > 0) {
                vc = sites.get(i);
            }
            results.add(new IBDLocus(seqs.get(i), starts.get(i), vc, states[i],
                    modelResults.statePosteriors[i], parentalIndicators.get(i),
                    observations.get(i).agreesWith(states[i]), observations.get(i)));
        }
        observations.clear();

        seqs.clear();
        starts.clear();
        sites.clear();
        parentalIndicators.clear();

        return results.iterator();
    }

    protected ModelResults runModel() {
        final int numObservations = observations.size();
        if (numObservations == 0) {
            final ModelResults nullResult = new ModelResults();
            nullResult.viterbiStates = new IBDState[0];
            nullResult.statePosteriors = new double[0][0];
            return nullResult;
        }

        final double[][] viterbi = new double[numObservations][IBDState.values().length];

        // forward, backward, and posterior probabilities for forward-backward algorithm
        final double[][] alphas = new double[numObservations][IBDState.values().length];
        final double[][] betas = new double[numObservations][IBDState.values().length];
        final double[][] gammas = new double[numObservations][IBDState.values().length];

        final short[][] path = new short[numObservations][IBDState.values().length];

        // forward and viterbi initialization
        for (final IBDState i : IBDState.values()) {
            viterbi[0][i.ordinal()] = initialProbs.get(i) + getEmission(i, observations.get(0));
            alphas[0][i.ordinal()] = initialProbs.get(i) + getEmission(i, observations.get(0));
        }

        // forward and viterbi
        for (int t = 1; t < numObservations; t++) {
            for (final IBDState j : IBDState.values()) {

                double maxProb = Double.NEGATIVE_INFINITY;
                int argMaxProb = -1;

                double sumProb = Double.NEGATIVE_INFINITY;

                for (final IBDState i : IBDState.values()) {
                    final double prob = viterbi[t - 1][i.ordinal()] + transitionProbs.get(i)[j.ordinal()] + getEmission(j, observations.get(t));
                    if (prob > maxProb) {
                        maxProb = prob;
                        argMaxProb = i.ordinal();
                    }

                    sumProb = MathUtils.approximateLog10SumLog10(sumProb, alphas[t - 1][i.ordinal()] + transitionProbs.get(i)[j.ordinal()] + getEmission(j, observations.get(t)));
                }
                viterbi[t][j.ordinal()] = maxProb;
                path[t][j.ordinal()] = (short) argMaxProb;

                alphas[t][j.ordinal()] = sumProb;
            }
        }

        // backward initialization
        for (final IBDState s : IBDState.values()) {
            betas[numObservations-1][s.ordinal()] = 0;
        }

        // backward
        for (int t = numObservations - 2; t >= 0; t--) {
            for (final IBDState i : IBDState.values()) {
                double sumProb = Double.NEGATIVE_INFINITY;
                for (final IBDState j : IBDState.values()) {
                    sumProb = MathUtils.approximateLog10SumLog10(sumProb, betas[t + 1][j.ordinal()] + transitionProbs.get(i)[j.ordinal()] + getEmission(j, observations.get(t + 1)));
                }
                betas[t][i.ordinal()] = sumProb;
            }
        }

        // backtrack through viterbi to find most likely states
        final IBDState[] states = new IBDState[numObservations];
        double maxProb = Double.NEGATIVE_INFINITY;
        int argMaxProb = -1;
        for (final IBDState s : IBDState.values()) {
            if (viterbi[numObservations-1][s.ordinal()] > maxProb) {
                maxProb = viterbi[numObservations-1][s.ordinal()];
                argMaxProb = s.ordinal();
            }
        }
        states[numObservations - 1] = IBDState.values()[argMaxProb];
        for (int i = numObservations - 1; i >= 1; i--) {
            states[i-1] = IBDState.values()[path[i][states[i].ordinal()]];
        }

        // gammas give posteriors of being in each state
        for (int t = 0; t < numObservations; t++) {
            double totalProb = Double.NEGATIVE_INFINITY;
            for (final IBDState i : IBDState.values()) {
                final double gamma = alphas[t][i.ordinal()] + betas[t][i.ordinal()];
                totalProb = MathUtils.approximateLog10SumLog10(totalProb, gamma);
                gammas[t][i.ordinal()] = gamma;
            }

            // normalize
            for (final IBDState j : IBDState.values()) {
                gammas[t][j.ordinal()] = gammas[t][j.ordinal()] - totalProb;
            }
        }

        final ModelResults modelResults = new ModelResults();
        modelResults.viterbiStates = states;
        modelResults.statePosteriors = gammas;
        return modelResults;
    }


    /**
     * Emission probabilities are simplistically modeled also, with (almost) uniform probabilities
     * for each allowable IBD observation in a given state and for error observations. A small amount
     * of extra probability is given to ZERO_OR_TWO observations (which arise when the quartet is all Het)
     * in IBD1 regions, which seem to be the most frequent error mode, likely arising from segmental duplications
     * or high-copy number regions
     *
     * @param s
     * @param o
     * @return
     */
    private Double getEmission(final IBDState s, final IBDObservation o) {
        switch (s) {
            case ZERO:
                if (o == IBDObservation.ZERO || o == IBDObservation.ZERO_OR_ONE || o == IBDObservation.ZERO_OR_TWO) {
                    return Math.log10(NON_ERROR_EMISSION_PROBABILITY / 3);
                } else {
                    return Math.log10(ERROR_EMISSION_PROBABILITY / 3);
                }
            case ONE:
                if (o == IBDObservation.ONE || o == IBDObservation.ZERO_OR_ONE || o == IBDObservation.ONE_OR_TWO) {
                    return Math.log10(NON_ERROR_EMISSION_PROBABILITY / 3);
                } else if (o == IBDObservation.ZERO_OR_TWO) {
                    return Math.log10(IBD1_HET_ERROR_EMISSION_PROBABILITY);
                } else {
                    return Math.log10(IDB1_OTHER_ERROR_EMISSION_PROBABILITY / 2);
                }
            case TWO:
                if (o == IBDObservation.TWO || o == IBDObservation.ZERO_OR_TWO || o == IBDObservation.ONE_OR_TWO) {
                    return Math.log10(NON_ERROR_EMISSION_PROBABILITY / 3);
                } else {
                    return Math.log10(ERROR_EMISSION_PROBABILITY / 3);
                }
        }
        return Double.NEGATIVE_INFINITY;
    }


    protected class ModelResults {
        IBDState[] viterbiStates;
        double[][] statePosteriors;
    }
}
