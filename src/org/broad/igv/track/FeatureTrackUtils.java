/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.track;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.tools.FeatureSearcher;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.CancellableProgressDialog;
import org.broad.igv.ui.util.IndefiniteProgressMonitor;
import org.broad.igv.ui.util.ProgressMonitor;
import org.broad.igv.util.LongRunningTask;
import org.broad.tribble.Feature;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.IOException;
import java.util.Iterator;

/**
 *
 */
class FeatureTrackUtils {

    /**
     * Find the next/previous feature which lies outside chr:initStart-initEnd
     * @param source
     * @param chr
     * @param initStart
     * @param initEnd
     * @param forward
     * @return
     * @throws IOException
     */
    public static Feature nextFeature(FeatureSource source, String chr, int initStart, int initEnd, boolean forward) throws IOException {
        Feature f = null;
        int binSize = source.getFeatureWindowSize();

        final Genome genome = GenomeManager.getInstance().getCurrentGenome();
        if (forward) {
            // Forward
            int nextStart = initEnd;
            String nextChr = chr;
            while (nextChr != null) {
                int chrLength = genome.getChromosome(nextChr).getLength();
                while (nextStart < chrLength) {
                    int nextEnd = binSize > 0 ? nextStart + source.getFeatureWindowSize() : chrLength;
                    Iterator<Feature> iter = source.getFeatures(nextChr, nextStart, nextEnd);
                    if (iter != null) {
                        // The check on position should not be necessary, but not all implementations of getFeatures
                        // obey the contract to return features only in the interval.
                        while (iter.hasNext()) {
                            Feature feat = iter.next();
                            if (feat.getStart() > nextStart) {
                                return feat;
                            }
                        }
                    }
                    nextStart = nextEnd;
                }
                nextChr = genome.getNextChrName(nextChr);
                nextStart = 0;
            }
        } else {
            // Reverse
            int nextEnd = initStart;
            String nextChr = chr;
            while (nextChr != null) {
                while (nextEnd > 0) {
                    int nextStart = binSize > 0 ? Math.max(0, nextEnd - source.getFeatureWindowSize()) : 0;
                    Iterator<Feature> iter = source.getFeatures(nextChr, nextStart, nextEnd);
                    if (iter != null && iter.hasNext()) {
                        // The check on position should not be necessary, but not all implementations of getFeatures
                        // obey the contract to return features only in the interval.
                        Feature prevFeature = null;
                        while (iter.hasNext()) {
                            Feature feat = iter.next();
                            if (feat.getStart() < nextEnd) {
                                prevFeature = feat;
                            }
                        }
                        if (prevFeature != null) {
                            return prevFeature;
                        }
                    }
                    nextEnd = nextStart;
                }
                nextChr = genome.getPrevChrName(nextChr);
                if (nextChr != null) {
                    nextEnd = genome.getChromosome(nextChr).getLength();
                }
            }
        }

        return f;
    }

    /**
     * Find the next/previous feature, using FeatureSearcher. A cancellable dialog will be displayed so the user can cancel the
     * search, null will be returned in that case.
     * @param source
     * @param chr
     * @param initStart
     * @param initEnd
     * @param forward
     * @param foundHandler Callback to invoke when searching is complete. Does not get called if search returns no result
     * @return
     * @throws IOException
     */
    public static void nextFeatureSearch(FeatureSource source, String chr, int initStart, int initEnd, boolean forward,
                                         final FeatureSearcher.IFeatureFound foundHandler) throws IOException{

        //Only allow one to be shown at a time
        if(CancellableProgressDialog.hasCancellableProgressDialog()){
            return;
        }

        Genome genome = GenomeManager.getInstance().getCurrentGenome();

        //We search backwards by setting a negative searchIncrement
        int searchIncrement = (initEnd - initStart) * (forward ? +1 : -1);
        int start = initStart + searchIncrement;

        IndefiniteProgressMonitor monitor = new IndefiniteProgressMonitor();

        final FeatureSearcher searcher = new FeatureSearcher(source, genome, chr, start, monitor);
        searcher.setSearchIncrement(searchIncrement);

        final ActionListener cancelListener = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                searcher.cancel();
            }
        };

        final CancellableProgressDialog dialog = CancellableProgressDialog.showCancellableProgressDialog(IGV.getMainFrame(), "Searching...", cancelListener, monitor);
        dialog.getProgressBar().setIndeterminate(true);

        monitor.addPropertyChangeListener(new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                if(evt.getPropertyName().equals(ProgressMonitor.PROGRESS_PROPERTY) &&  (Integer) evt.getNewValue() >= 100){
                    Iterator<? extends Feature> result = searcher.getResult();
                    if(result != null) foundHandler.processResult(result);
                }
            }
        });

        LongRunningTask.submit(searcher);
    }
}
