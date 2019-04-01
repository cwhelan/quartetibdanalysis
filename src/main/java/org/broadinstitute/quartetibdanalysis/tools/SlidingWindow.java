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
import java.util.HashMap;
import java.util.Map;

/**
 * Created by cwhelan on 5/1/14.
 *
 * A sliding window that keeps counts of objects of type K in that window
 *
 */
public class SlidingWindow<K> {

    private int windowSize;

    private Deque<K> window;
    private Deque<Integer> coords;
    private final Map<K, Integer> windowCounts = new HashMap<>();

    public SlidingWindow(final int windowSize) {
        this.windowSize = windowSize;
        window = new ArrayDeque<>(windowSize);
        coords = new ArrayDeque<>(windowSize / 2);
    }

    /**
     * Add an observation of type K at coordinate start
     * @param start
     * @param key
     * @return the counts of the center of the current window or null if we don't have enough observations to complete
     * a window yet
     */
    public WindowCenterCount next(final int start, final K key) {
        incrementCount(start, key);
        int windowCenterCoord = 0;
        if (window.size() > windowSize / 2) {
            windowCenterCoord = coords.removeLast();
        }

        if (window.size() < windowSize) {
            return null;
        }

        final Map<K, Double> returnCounts = new HashMap<>();
        for (final K c : windowCounts.keySet()) {
            returnCounts.put(c, ((double) windowCounts.get(c)) / windowSize);
        }
        decrementCount();
        final WindowCenterCount result = new WindowCenterCount();
        result.start = windowCenterCoord;
        result.counts = returnCounts;
        return result;
    }


    private void incrementCount(final int start, final K key) {
        window.addFirst(key);
        coords.addFirst(start);
        if (! windowCounts.containsKey(key)) {
            windowCounts.put(key, 1);
        } else {
            windowCounts.put(key, windowCounts.get(key) + 1);
        }
    }

    private void decrementCount() {
        final K last = window.removeLast();
        windowCounts.put(last, windowCounts.get(last) - 1);
    }

    public class WindowCenterCount {
        int start;
        private Map<K, Double> counts;

        public Double getCount(final K key) {
            if (counts.containsKey(key)) {
                return counts.get(key);
            } else {
                return 0.0;
            }
        }
    }
}
