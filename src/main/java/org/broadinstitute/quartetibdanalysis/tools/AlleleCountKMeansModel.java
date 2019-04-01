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

import java.util.HashMap;
import java.util.Map;

/**
* Created by cwhelan on 5/6/14.
 *
 * Using a k-means model, this class determines the most likely IBD state
 * at a locus given the proportion of shared minor alleles in a window
 * centered on that locus.
 *
*/
public class AlleleCountKMeansModel {
    private Map<SharedMinorAlleleClass, Double> means = new HashMap<>();
    private Map<SharedMinorAlleleClass, Double> sds = new HashMap<>();
    private Map<SharedMinorAlleleClass, Double> ibd0ClusterMean = new HashMap<>();
    private Map<SharedMinorAlleleClass, Double> ibd1ClusterMean = new HashMap<>();
    private Map<SharedMinorAlleleClass, Double> ibd2ClusterMean = new HashMap<>();

    public AlleleCountKMeansModel() {
        // These were from a K-means model trained on a human data set using results from ISCA
        // as training labels. The model was built in R using the kmeans() function from
        // the stats package. This model may not generalize well to other data sets, and therefore
        // should be considered experimental and not ready for use on production data
        // without further testing and verification.
        means.put(SharedMinorAlleleClass.HOMVAR_HOMVAR, 0.05033174);
        means.put(SharedMinorAlleleClass.HETVAR_HETVAR, 0.36483247);
        means.put(SharedMinorAlleleClass.HOMVAR_HETVAR, 0.12042467);
        means.put(SharedMinorAlleleClass.HETVAR_HOMREF, 0.44547881);
        means.put(SharedMinorAlleleClass.HOMVAR_HOMREF, 0.01943231);

        sds.put(SharedMinorAlleleClass.HOMVAR_HOMVAR, 0.06060862);
        sds.put(SharedMinorAlleleClass.HETVAR_HETVAR, 0.19861271);
        sds.put(SharedMinorAlleleClass.HOMVAR_HETVAR, 0.07279400);
        sds.put(SharedMinorAlleleClass.HETVAR_HOMREF, 0.17219448);
        sds.put(SharedMinorAlleleClass.HOMVAR_HOMREF, 0.03922841);

        ibd0ClusterMean.put(SharedMinorAlleleClass.HOMVAR_HOMVAR, -0.60162856);
        ibd0ClusterMean.put(SharedMinorAlleleClass.HETVAR_HETVAR, -0.8137241);
        ibd0ClusterMean.put(SharedMinorAlleleClass.HOMVAR_HETVAR, -0.2495874);
        ibd0ClusterMean.put(SharedMinorAlleleClass.HETVAR_HOMREF, 0.90989689);
        ibd0ClusterMean.put(SharedMinorAlleleClass.HOMVAR_HOMREF, 1.5185183);

        ibd1ClusterMean.put(SharedMinorAlleleClass.HOMVAR_HOMVAR, 0.04409738);
        ibd1ClusterMean.put(SharedMinorAlleleClass.HETVAR_HETVAR, -0.1544063);
        ibd1ClusterMean.put(SharedMinorAlleleClass.HOMVAR_HETVAR, 0.4203370);
        ibd1ClusterMean.put(SharedMinorAlleleClass.HETVAR_HOMREF, 0.09125241);
        ibd1ClusterMean.put(SharedMinorAlleleClass.HOMVAR_HOMREF, -0.4669271);

        ibd2ClusterMean.put(SharedMinorAlleleClass.HOMVAR_HOMVAR, 0.83843249);
        ibd2ClusterMean.put(SharedMinorAlleleClass.HETVAR_HETVAR, 2.1130272);
        ibd2ClusterMean.put(SharedMinorAlleleClass.HOMVAR_HETVAR, -1.4910337);
        ibd2ClusterMean.put(SharedMinorAlleleClass.HETVAR_HOMREF, -1.99044294);
        ibd2ClusterMean.put(SharedMinorAlleleClass.HOMVAR_HOMREF, -0.4896625);
    }

    public Map<SharedMinorAlleleClass, Double> scale(final SlidingWindow<SharedMinorAlleleClass>.WindowCenterCount counts) {
        final Map<SharedMinorAlleleClass, Double> results = new HashMap<>();
        for (SharedMinorAlleleClass c : means.keySet()) {
            results.put(c, (counts.getCount(c) - means.get(c)) / sds.get(c));
        }
        return results;
    }

    public IBDState ibdClass(final SlidingWindow<SharedMinorAlleleClass>.WindowCenterCount counts) {
        final Map<SharedMinorAlleleClass, Double> scaled = scale(counts);
        final double dist0 = dist(scaled, ibd0ClusterMean);
        final double dist1 = dist(scaled, ibd1ClusterMean);
        final double dist2 = dist(scaled, ibd2ClusterMean);
        if (dist0 < dist1 && dist0 < dist2) {
            return IBDState.ZERO;
        }
        if (dist1 < dist0 && dist1 < dist2) {
            return IBDState.ONE;
        }
        return IBDState.TWO;
    }

    private double dist(final Map<SharedMinorAlleleClass, Double> x, final Map<SharedMinorAlleleClass, Double> clusterMean) {
        double distance = 0;
        for (SharedMinorAlleleClass c : clusterMean.keySet()) {
            distance += Math.pow(x.get(c) - clusterMean.get(c), 2);
        }
        return Math.sqrt(distance);
    }
}
