package org.broad.igv.util.liftover;

import htsjdk.tribble.Feature;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.util.Interval;
import org.broad.igv.util.IntervalTree;

import java.util.ArrayList;
import java.util.List;

public class Chain {
    String tName;
    int tSize;
    int tStart;
    int tEnd;
    String qName;
    int qSize;
    int qStart;
    int qEnd;
    String id;

    IntervalTree<int[]> tree;

    public Chain(String tName, int tSize, int tStart, int tEnd, String qName, int qSize, int qStart, int qEnd, String id) {
        this.tName = tName;
        this.tSize = tSize;
        this.tStart = tStart;
        this.tEnd = tEnd;
        this.qName = qName;
        this.qSize = qSize;
        this.qStart = qStart;
        this.qEnd = qEnd;
        this.id = id;
    }


    public void setAlignments(List<String[]> lines) {

        this.tree = new IntervalTree();

        int tStart = this.tStart;
        int qStart = this.qStart;

        for (String[] line : lines) {

            int size = Integer.parseInt(line[0]);
            int tEnd = tStart + size;
            int qEnd = qStart + size;

            int[] value = {qStart, qEnd};
            this.tree.insert(new Interval<>(tStart, tEnd, value));

            if(line.length == 3) {
                tStart += size + Integer.parseInt(line[1]);
                qStart += size + Integer.parseInt(line[2]);
            }

        }

    }

    /**
     * Map a span represented by the array [start, end] in target coordinates to query coordinates
     *
     * @param span
     * @return
     */
    public List<int[]> map(int[] span) {

        List<int[]> mapped = new ArrayList<>();
        List<Interval<int[]>> intervals = this.tree.findOverlapping(span[0], span[1]);
        if(intervals != null) {
            for(Interval<int[]> interval : intervals) {

                int ds = span[0] - interval.getLow();

                int [] qspan = interval.getValue();
                int start = Math.max(qspan[0], qspan[0] + ds);
                int end = Math.min(qspan[1], qspan[0] + (span[1] - span[0]) + ds);
                int [] mappedSpan = {start, end};
                mapped.add(mappedSpan);
            }
        }

        return mapped;
    }
}
