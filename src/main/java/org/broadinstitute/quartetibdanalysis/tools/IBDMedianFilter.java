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

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.EnumMap;
import java.util.Map;

/**
 * Created by cwhelan on 5/6/14.
 *
 * Given a window size, applies a median filter, setting the value of the
 * locus at the center of the window to the median of the values in the
 * window.
 *
 */
public class IBDMedianFilter {
    private int windowSize;

    private Deque<IBDState> ibdClassWindow;
    private Deque<Integer> coords;

    private final Map<IBDState, Integer> clusterCounts = new EnumMap<>(IBDState.class);

    public IBDMedianFilter(final int windowSize) {
        this.windowSize = windowSize;
        ibdClassWindow = new ArrayDeque<>(windowSize);
        coords = new ArrayDeque<>(windowSize / 2);
        clusterCounts.put(IBDState.ZERO, 0);
        clusterCounts.put(IBDState.ONE, 0);
        clusterCounts.put(IBDState.TWO, 0);
    }

    public FilteredValue next(final int start, final IBDState ibdClass) {
        incrementCount(start, ibdClass);
        int windowCenterCoord = 0;
        if (ibdClassWindow.size() > windowSize / 2) {
            windowCenterCoord = coords.removeLast();
        }

        if (ibdClassWindow.size() < windowSize) {
            return null;
        }

        final IBDState median = median();
        decrementCount();
        final FilteredValue result = new FilteredValue();
        result.start = windowCenterCoord;
        result.ibdClass = median;
        return result;
    }

    private IBDState median() {
        if (clusterCounts.get(IBDState.ZERO) >= windowSize / 2) {
            return IBDState.ZERO;
        } else if (clusterCounts.get(IBDState.ZERO) + clusterCounts.get(IBDState.ONE) >= windowSize / 2) {
            return IBDState.ONE;
        } else {
            return IBDState.TWO;
        }
    }


    private void incrementCount(final int start, final IBDState next) {
        ibdClassWindow.addFirst(next);
        coords.addFirst(start);
        clusterCounts.put(next, clusterCounts.get(next) + 1);
    }

    private void decrementCount() {
        final IBDState last = ibdClassWindow.removeLast();
        clusterCounts.put(last, clusterCounts.get(last) - 1);
    }

    public class FilteredValue {
        int start;
        IBDState ibdClass;
    }

}
