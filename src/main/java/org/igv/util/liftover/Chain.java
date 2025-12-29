package org.igv.util.liftover;

import htsjdk.tribble.Feature;
import org.igv.feature.IGVFeature;
import org.igv.feature.Range;
import org.igv.util.Interval;
import org.igv.util.IntervalTree;

import java.util.ArrayList;
import java.util.List;

/**
 * Represents a single "chain" from a UCSC chain format file (see https://genome.ucsc.edu/goldenPath/help/chain.html)
 */
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


    /**
     * Set the pairwise alignments from alignment lines of the chain file.
     * @param lines
     */
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
     * Map a region in target coordinates to query coordinates
     *
     * @param span
     * @return
     */
    public List<Range> map(Range span) {

        List<Range> mapped = new ArrayList<>();

        if(span.getChr().equals(this.tName)) {
            List<Interval<int[]>> intervals = this.tree.findOverlapping(span.start, span.end);
            if (intervals != null) {
                for (Interval<int[]> interval : intervals) {

                    int ds = span.start - interval.getLow();

                    int[] qspan = interval.getValue();
                    int start = Math.max(qspan[0], qspan[0] + ds);
                    int end = Math.min(qspan[1], qspan[0] + (span.getLength()) + ds);
                    Range mappedSpan = new Range(this.qName, start, end);
                    mapped.add(mappedSpan);
                }
            }
        }

        return mapped;
    }
}
