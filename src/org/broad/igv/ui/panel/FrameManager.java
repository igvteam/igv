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

package org.broad.igv.ui.panel;

import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Locus;
import org.broad.igv.lists.GeneList;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.action.SearchCommand;
import org.broad.igv.ui.util.MessageUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * @author jrobinso
 * @date Sep 10, 2010
 */
public class FrameManager {

    private static List<ReferenceFrame> frames = new ArrayList();
    private static ReferenceFrame defaultFrame;

    public static final String DEFAULT_FRAME_NAME = "genome";

    static {
        frames.add(getDefaultFrame());
    }

    public synchronized static ReferenceFrame getDefaultFrame() {
        if (defaultFrame == null) {
            defaultFrame = new ReferenceFrame(DEFAULT_FRAME_NAME);
        }
        return defaultFrame;
    }

    public static List<ReferenceFrame> getFrames() {
        return frames;
    }

    public static ReferenceFrame getFrame(String frameName){
        for(ReferenceFrame frame: frames){
            if(frame.getName().equals(frameName)){
                return frame;
            }
        }
        return null;
    }

    public static void setFrames(List<ReferenceFrame> f) {
        frames = f;
    }

    public static boolean isGeneListMode() {
        return frames.size() > 1;
    }

    public static int getStateHash() {
        if(isGeneListMode()) {
            String hs = "";
            for(ReferenceFrame frame : frames) {
                hs = hs + frame.getStateHash();
            }
            return hs.hashCode();
        }
        else {
            return defaultFrame.getStateHash();
        }
    }

    public static void setToDefaultFrame(String searchString) {
        frames.clear();
        if (searchString != null) {
            Locus locus = getLocus(searchString, 0);
            if (locus != null) {
                getDefaultFrame().jumpTo(locus);
            }
        }
        frames.add(getDefaultFrame());
        getDefaultFrame().recordHistory();
    }

    private static boolean addNewFrame(String searchString){
        boolean locusAdded = false;
        Locus locus = getLocus(searchString);
        if (locus != null) {
            ReferenceFrame referenceFrame = new ReferenceFrame(searchString);
            referenceFrame.jumpTo(locus);
            locusAdded = frames.add(referenceFrame);
        }
        return locusAdded;
    }

    public static void resetFrames(GeneList gl) {

        frames.clear();

        if (gl == null) {
            frames.add(getDefaultFrame());
        } else {
            List<String> lociNotFound = new ArrayList();
            List<String> loci = gl.getLoci();
            if (loci.size() == 1) {
                Locus locus = getLocus(loci.get(0));
                if (locus == null) {
                    lociNotFound.add(loci.get(0));
                } else {
                    IGV.getInstance().getSession().setCurrentGeneList(null);
                    getDefaultFrame().jumpTo(locus.getChr(), locus.getStart(), locus.getEnd());
                }
            } else {
                for (String searchString : gl.getLoci()) {
                    if(!addNewFrame(searchString)){
                        lociNotFound.add(searchString);
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


    /**
     * Uses default flanking region with
     * {@link #getLocus(String, int)}
     * @param searchString
     * @return
     */
    public static Locus getLocus(String searchString) {
        int flankingRegion = PreferenceManager.getInstance().getAsInt(PreferenceManager.FLANKING_REGION);
        return getLocus(searchString, flankingRegion);
    }

    /**
     * Runs a search for the specified string, and returns a locus
     * of the given region with additional space on each side.
     * Note: We DO NOT add the flanking region if the {@code searchString}
     * is a locus (e.g. chr1:50-100), only if it's a gene or feature name (or something else)
     *
     * @param searchString
     * @param flankingRegion
     * @return The found locus, null if not found
     */
    public static Locus getLocus(String searchString, int flankingRegion) {
        SearchCommand cmd = new SearchCommand(getDefaultFrame(), searchString);
        List<SearchCommand.SearchResult> results = cmd.runSearch(searchString);
        Locus locus = null;
        for (SearchCommand.SearchResult result : results) {
            if (result.getType() != SearchCommand.ResultType.ERROR) {
                int delta = 0;

                if (result.getType() != SearchCommand.ResultType.LOCUS) {
                    if (flankingRegion < 0) {
                        delta = (-flankingRegion * (result.getEnd() - result.getStart())) / 100;
                    } else {
                        delta = flankingRegion;
                    }
                }

                int start = result.getStart() - delta;
                //Don't allow flanking region to extend past origin
                //There are some circumstances in which we render before origin (e.g. soft-clips)
                //so we are conservative
                if (start < 0 && result.getStart() >= -1) {
                    start = 0;
                }
                locus = new Locus(
                        result.getChr(),
                        start,
                        result.getEnd() + delta);
                //We just take the first result
                break;
            }
        }
        return locus;
    }

    public static void removeFrame(ReferenceFrame frame) {
        frames.remove(frame);
    }

    public static void sortFrames(final Track t) {

        Collections.sort(frames, new Comparator<ReferenceFrame>() {
            @Override
            public int compare(ReferenceFrame o1, ReferenceFrame o2) {
                float s1 = t.getRegionScore(o1.getChromosome().getName(), (int) o1.getOrigin(), (int) o1.getEnd(),
                        o1.getZoom(), RegionScoreType.SCORE, o1.getName());
                float s2 = t.getRegionScore(o2.getChromosome().getName(), (int) o2.getOrigin(), (int) o2.getEnd(),
                        o2.getZoom(), RegionScoreType.SCORE, o2.getName());
                return (s1 == s2 ? 0 : (s1 > s2) ? -1 : 1);
            }
        });

    }

    /**
     * Increment the zoom level of the visibile frame(s). Supports batch commands zoomIn and zoomOut
     *
     * @param zoom the zoom level increment, usually -1 or 1
     */
    public static void incrementZoom(int zoom) {

        if(isGeneListMode()) {
            for(ReferenceFrame frame : getFrames()) {
                frame.doZoomIncrement(zoom);
            }
        }
        else {
            getDefaultFrame().doZoomIncrement(zoom);
        }
    }


}

