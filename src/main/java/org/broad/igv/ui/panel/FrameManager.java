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

import org.broad.igv.event.GenomeChangeEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.genome.ChromosomeNameComparator;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.lists.GeneList;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.session.Session;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.action.SearchCommand;
import org.broad.igv.ui.util.MessageUtils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * @author jrobinso
 * @date Sep 10, 2010
 */
public class FrameManager implements IGVEventObserver {

    private static List<ReferenceFrame> frames = new ArrayList();
    private static ReferenceFrame defaultFrame;

    public static final String DEFAULT_FRAME_NAME = "genome";

    static {
        frames.add(getDefaultFrame());
    }

    public static String getCurrentLocusString() {
        if (isGeneListMode()) {
            String name = frames.get(0).getName();
            for (int i = 1; i < frames.size(); i++) {
                name += " | " + frames.get(i).getName();
            }
            return name;
        } else {
            return defaultFrame.getFormattedLocusString();
        }
    }

    public synchronized static ReferenceFrame getDefaultFrame() {
        if (defaultFrame == null) {
            defaultFrame = new ReferenceFrame(DEFAULT_FRAME_NAME);
        }
        return defaultFrame;
    }


    public static ReferenceFrame getFirstFrame() {
        return isGeneListMode() ? frames.get(0) : defaultFrame;
    }

    public static List<ReferenceFrame> getFrames() {
        return frames;
    }

    public static ReferenceFrame getFrame(String frameName) {
        for (ReferenceFrame frame : frames) {
            if (frame.getName().equals(frameName)) {
                return frame;
            }
        }
        return null;
    }

    public static void setFrames(List<ReferenceFrame> f) {
        frames = f;
        IGVEventBus.getInstance().post(new ChangeEvent(frames));
    }

    public static boolean isGeneListMode() {
        return frames.size() > 1;
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

        IGVEventBus.getInstance().post(new ChangeEvent(frames));
    }

    private static boolean addNewFrame(String searchString) {
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
                    frames.add(getDefaultFrame());
                }
            } else {
                for (String searchString : gl.getLoci()) {
                    if (!addNewFrame(searchString)) {
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

        IGVEventBus.getInstance().post(new ChangeEvent(frames));
    }

    /**
     * @return The minimum scale among all active frames
     * TODO -- track this with "rescale" events, rather than compute on the fly
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
     *
     * @param searchString
     * @return
     */
    public static Locus getLocus(String searchString) {
        int flankingRegion = PreferencesManager.getPreferences().getAsInt(Constants.FLANKING_REGION);
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

        if (isGeneListMode()) {
            for (ReferenceFrame frame : getFrames()) {
                frame.doZoomIncrement(zoom);
            }
        } else {
            getDefaultFrame().doZoomIncrement(zoom);
        }
    }

    /**
     * Add a new list of loci (frames).  First, a check is made for overlap with current frames, any any loci overlap
     * the current frame is expanded.  If no overlaps are found a new frame is created.
     * <p>
     * This method was added so support "circular view" interactions.
     *
     * @param newLociStrings
     */
    public static void addFrames(List<String> newLociStrings) {

        if (newLociStrings.size() != 2) {
            throw new RuntimeException("Unexpected list size: " + newLociStrings.size());
        }

        // Convert loci strings to locus objects
        List<Locus> newLoci = newLociStrings.stream()
                .map(s -> {
                    Locus locus = Locus.fromString(s);
                    locus.chr = GenomeManager.getInstance().getCurrentGenome().getCanonicalChrName(locus.chr);
                    return locus;
                })
                .collect(Collectors.toList());

        // If new loci overlap combine
        if(newLoci.get(0).overlaps(newLoci.get(1))) {
            newLoci.get(0).start = Math.min(newLoci.get(0).start, newLoci.get(1).start);
            newLoci.get(0).end = Math.max(newLoci.get(0).end, newLoci.get(1).end);
            newLoci.remove(1);
        }

        final List<ReferenceFrame> currentFrames =
                isGeneListMode() ? frames :
                        defaultFrame.chrName.equals("All") ?
                                Arrays.asList() :
                                Arrays.asList(defaultFrame);


        Set<ReferenceFrame> usedFrames = new HashSet<>();

        List<String> loci = new ArrayList<>();
        for (Locus locus : newLoci) {
            boolean found = false;
            for (ReferenceFrame ref : currentFrames) {
                if (locus != null && ref.getChrName().equals(locus.chr) && ref.getCurrentRange().overlaps(locus)) {
                    Range union = ref.getCurrentRange().union(locus);
                    loci.add(Locus.getFormattedLocusString(union.getChr(), union.getStart(), union.getEnd()));
                    found = true;
                    usedFrames.add(ref);
                    break;
                }
            }
            if (!found) {
                loci.add(Locus.getFormattedLocusString(locus.getChr(), locus.getStart(), locus.getEnd()));
            }
        }

        for (ReferenceFrame ref : currentFrames) {
            if (!usedFrames.contains(ref)) {
                loci.add(ref.getFormattedLocusString());
            }
        }

        GeneList geneList = new GeneList("Current frames", loci, false);
        Session currentSession = IGV.getInstance().getSession();
        currentSession.setCurrentGeneList(geneList);

        // sort the frames by position
        currentSession.sortGeneList((n0, n1) -> {
            ReferenceFrame f0 = getFrame(n0);
            ReferenceFrame f1 = getFrame(n1);

            String chr0 = f0 == null ? "" : f0.getChrName();
            String chr1 = f1 == null ? "" : f1.getChrName();
            int s0 = f0 == null ? 0 : f0.getCurrentRange().getStart();
            int s1 = f1 == null ? 0 : f1.getCurrentRange().getStart();

            int chrComp = ChromosomeNameComparator.get().compare(chr0, chr1);
            if (chrComp != 0) return chrComp;
            return s0 - s1;
        });
        IGV.getInstance().resetFrames();
    }


    @Override
    public void receiveEvent(Object event) {
        if (event instanceof GenomeChangeEvent) {
            Genome newGenome = ((GenomeChangeEvent) event).genome;
            boolean force = true;
            getDefaultFrame().setChromosomeName(newGenome.getHomeChromosome(), force);
        }
    }


    public static class ChangeEvent {
        List<ReferenceFrame> frames;

        public ChangeEvent(List<ReferenceFrame> frames) {
            this.frames = frames;
        }

        public List<ReferenceFrame> getFrames() {
            return frames;
        }
    }

}

