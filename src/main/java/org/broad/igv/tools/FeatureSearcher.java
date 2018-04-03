/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.tools;

import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.IndefiniteProgressMonitor;
import htsjdk.tribble.Feature;

import java.io.IOException;
import java.util.Iterator;

/**
 * Used for searching for the next feature, given a source.
 * We simply call getFeature on that source repeatedly over stepped windows
 * until we find something
 * User: jacob
 * Date: 2013-Feb-21
 */
public class FeatureSearcher implements Runnable {

    private static Logger log = Logger.getLogger(FeatureSearcher.class);

    private FeatureTrack track = null;
    private FeatureSource<? extends Feature> source = null;

    private static final int DEFAULT_SEARCH_INCREMENT = 100000;

    /**
     * After searching a window, the increment over which to move.
     * If smaller than searchWindowSize, it will be a sliding search with overlap.
     * If larger, there will be gaps. Make negative to search backwards
     */
    private int searchIncrement = DEFAULT_SEARCH_INCREMENT;

    /**
     * The window size over which to search, in base pairs
     */
    private int searchWindowSize = searchIncrement;

    private volatile Iterator<? extends Feature> result = null;

    private volatile boolean isRunning = false;
    private volatile boolean wasCancelled = false;
    private volatile boolean done = false;

    private final Genome genome;
    private String chr;
    private int start;
    private int end;
    private IndefiniteProgressMonitor monitor;

    public FeatureSearcher(FeatureSource<? extends Feature> source, Genome genome, String chr, int start){
        this(source, genome, chr, start, null);
    }

    /**
     *
     * @param source FeatureSource which we are searching
     * @param genome
     * @param chr
     * @param start
     * @param monitor Optional (may be null)
     */
    public FeatureSearcher(FeatureSource<? extends Feature> source, Genome genome, String chr, int start, IndefiniteProgressMonitor monitor){
        assert source != null;
        this.source = source;
        this.genome = genome;
        this.monitor = monitor;
        this.initSearchCoords(chr, start);
    }

    private void initSearchCoords(String chr, int start){
        this.chr = chr;
        this.start = start;
        this.end = start + searchWindowSize;
    }

    private void incrementSearchCoords(){
        this.start += searchIncrement;
        int maxCoord = Integer.MAX_VALUE - searchWindowSize;
        int minCoord = 0;

        if(this.genome != null){
            maxCoord = genome.getChromosome(chr).getLength();
        }

        boolean outsideBounds = start >= maxCoord || start < minCoord;


        if (outsideBounds) {
            String lastChr = chr;
            chr = null;
            if(genome != null){
                if(start >= maxCoord){
                    chr = genome.getNextChrName(lastChr);
                }else if(start < minCoord){
                    chr = genome.getPrevChrName(lastChr);
                }
            }

            if (chr == null) {
                //No next chromosome, done searching
                start = end = -1;
                this.cancel();
                return;
            } else {
                maxCoord = genome.getChromosome(chr).getLength();
                start = searchIncrement > 0 ? minCoord : maxCoord - searchWindowSize;
            }
        }
        this.end = start + searchWindowSize;
        this.end = Math.min(this.end, maxCoord);
    }

    private Iterator<? extends Feature> getFeatures(String chr, int start, int end) throws IOException{
        if(track != null) return track.getFeatures(chr, start, end).iterator();
        if(source != null) return source.getFeatures(chr, start, end);
        throw new IllegalStateException("Have no FeatureTrack or FeatureSource from which to get features");
    }

    /**
     * Signal the searcher to stop. Note that stopping may not be instantaneous
     */
    public void cancel(){
        this.wasCancelled = true;
    }

    public boolean isDone(){
        return this.done;
    }

    public Iterator<? extends Feature> getResult(){
        if(this.isRunning) return null;
        return this.result;
    }

    /**
     * Set the search increment, can be either positive or negative (to search backwards).
     * Also sets the searchWindowSize
     * @param searchIncrement
     */
    public void setSearchIncrement(int searchIncrement) {
        if(this.isRunning) throw new IllegalStateException("Cannot set search increment while searching");
        this.searchIncrement = searchIncrement;
        this.searchWindowSize = Math.abs(searchIncrement);
        this.end = this.start + this.searchWindowSize;
    }

    @Override
    public void run() {
        isRunning = true;
        Iterator<? extends Feature> rslt = null;
        int counter = 0;
        int updateInterval = 20 * (int) (DEFAULT_SEARCH_INCREMENT / (1.0 * searchIncrement) );
        //Keep updateInterval above 0
        updateInterval = Math.max(updateInterval, 100);

        if(this.monitor != null){
            this.monitor.start();
        }

        while(isRunning && !wasCancelled){
            try {
                if(this.monitor != null){
                    if(counter == 0){
                        String status = String.format("Searching: %s:%d-%d", chr, start, end);
                        this.monitor.updateStatus(status);
                    }
                    counter = (counter + 1) % updateInterval;
                }
                rslt = getFeatures(chr, start, end);
                if(rslt != null && rslt.hasNext()){
                    //Found something
                    this.result = rslt;
                    break;
                }else{
                    //Didn't find anything, keep going
                    incrementSearchCoords();
                }
            } catch (IOException e) {
                log.error("Error searching for feature", e);
                break;
            }
        }
        this.isRunning = false;
        this.done = true;
        if(this.monitor != null){
            this.monitor.stop();
        }
    }

    /**
     * Listener for handling search result
     */
    public static interface IFeatureFound{
        void processResult(Iterator<? extends Feature> searchResult);
    }

    public static class GotoFeatureHandler implements IFeatureFound{
        @Override
        public void processResult(Iterator<? extends Feature> searchResult) {
            ReferenceFrame frame = FrameManager.getDefaultFrame();
            Feature f = searchResult.next();

            String chr = GenomeManager.getInstance().getCurrentGenome().getCanonicalChrName(f.getChr());
            double newCenter = f.getStart();
            if (!chr.equals(frame.getChrName())) {
                // Switch chromosomes.  We have to do some tricks to maintain the same resolution scale.
                double range = frame.getEnd() - frame.getOrigin();
                int newOrigin = (int) Math.max(newCenter - range / 2, 0);
                int newEnd = (int) (newOrigin + range);
                frame.jumpTo(chr, newOrigin, newEnd);
            } else {
                frame.centerOnLocation(newCenter);
            }
        }
    }
}
