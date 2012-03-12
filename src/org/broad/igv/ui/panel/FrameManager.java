/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.ui.panel;

import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.NamedFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.lists.GeneList;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.action.SearchCommand;
import org.broad.igv.ui.util.MessageUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 * @date Sep 10, 2010
 */
public class FrameManager {

    private static List<ReferenceFrame> frames = new ArrayList();
    private static ReferenceFrame defaultFrame;

    static {
        defaultFrame = new ReferenceFrame("genome");
        frames.add(defaultFrame);
    }

    public static ReferenceFrame getDefaultFrame() {
        return defaultFrame;
    }


    public static List<ReferenceFrame> getFrames() {
        return frames;
    }

    public static void setFrames(List<ReferenceFrame> f) {
        frames = f;
    }

    public static boolean isGeneListMode() {
        return frames.size() > 1;
    }


    public static void setToDefaultFrame(String searchString) {
        frames.clear();
        if (searchString != null) {
            Locus locus = getLocus(searchString, 0);
            if (locus != null) {
                defaultFrame.setInterval(locus);
            }
        }
        frames.add(defaultFrame);
        defaultFrame.recordHistory();
    }


    public static void resetFrames(GeneList gl) {

        frames.clear();

        if (gl == null) {
            frames.add(defaultFrame);
        } else {
            int flankingRegion = PreferenceManager.getInstance().getAsInt(PreferenceManager.FLANKING_REGION);
            List<String> lociNotFound = new ArrayList();
            List<String> loci = gl.getLoci();
            if (loci.size() == 1) {
                Locus locus = getLocus(loci.get(0), flankingRegion);
                if (locus == null) {
                    lociNotFound.add(loci.get(0));
                } else {
                    IGV.getInstance().getSession().setCurrentGeneList(null);
                    defaultFrame.jumpTo(locus.getChr(), locus.getStart(), locus.getEnd());
                }
            } else {
                for (String searchString : gl.getLoci()) {
                    Locus locus = getLocus(searchString, flankingRegion);
                    if (locus == null) {
                        lociNotFound.add(searchString);
                    } else {
                        ReferenceFrame referenceFrame = new ReferenceFrame(searchString);
                        referenceFrame.setInterval(locus);
                        frames.add(referenceFrame);
                    }
                }
            }

            if (lociNotFound.size() > 1) {
                StringBuffer message = new StringBuffer();
                message.append("<html>The following loci could not be found in the currently loaded annotation sets: <br>");
                for (String s : lociNotFound) {
                    message.append(s + " ");
                }
                MessageUtils.showMessage(message.toString());

            }
        }
    }

    /**
     * @return The minimum scale among all active frames
     *         TODO -- track this with "rescale" events, rather than compute on the fly
     */
    public static double getMinimumScale() {
        double minScale = Double.MAX_VALUE;
        for (ReferenceFrame frame : frames) {
            minScale = Math.min(minScale, frame.getScale());
        }
        return minScale;
    }


    public static Locus getLocus(String name) {
        int flankingRegion = PreferenceManager.getInstance().getAsInt(PreferenceManager.FLANKING_REGION);
        return getLocus(name, flankingRegion);
    }

    /**
     * This will actually result in a double call when used from search box,
     * because that calls to here as well. Not inifinite, won't overflow, just
     * a waste of time.
     *
     * @param searchString
     * @param flankingRegion
     * @return
     */
    public static Locus getLocusNew(String searchString, int flankingRegion) {
        SearchCommand cmd = new SearchCommand(getDefaultFrame(), searchString);
        List<SearchCommand.SearchResult> results = cmd.runSearch(searchString);
        Locus locus = null;
        for (SearchCommand.SearchResult result : results) {
            if (result.getType() != SearchCommand.ResultType.ERROR) {
                locus = new Locus(
                        result.getChr(),
                        result.getStart() - flankingRegion,
                        result.getEnd() + flankingRegion);
                //We just take the first result
                break;
            }
        }
        return locus;
    }


    /**
     * TODO This duplicates functionality in SearchCommand. Should merge these
     *
     * @param searchString
     * @param flankingRegion
     * @return
     */
    public static Locus getLocus(String searchString, int flankingRegion) {

        Locus locus = null;
        NamedFeature feature = FeatureDB.getFeature(searchString.toUpperCase().trim());
        if (feature != null) {
            locus = new Locus(
                    feature.getChr(),
                    feature.getStart() - flankingRegion,
                    feature.getEnd() + flankingRegion);
        } else {
            locus = new Locus(searchString);
            String chr = locus.getChr();
            if (chr != null) {
                return locus;
            } else {
                locus = null;
                if (IGV.hasInstance()) {
                    Genome genome = IGV.getInstance().getGenomeManager().getCurrentGenome();
                    if (genome != null) {
                        Chromosome chromsome = genome.getChromosome(searchString);
                        if (chromsome != null) {
                            locus = new Locus(chromsome.getName(), 0, chromsome.getLength());
                        }
                    }
                }
            }
        }
        return locus;
    }

    public static void removeFrame(ReferenceFrame frame) {
        frames.remove(frame);
    }


    public static void reset(String chr) {
        setToDefaultFrame(null);
        getDefaultFrame().setChrName(chr);
        getDefaultFrame().computeMaxZoom();
        getDefaultFrame().invalidateLocationScale();
    }

}

