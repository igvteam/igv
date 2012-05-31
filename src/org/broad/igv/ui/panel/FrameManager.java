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

package org.broad.igv.ui.panel;

import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.NamedFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.xome.Block;
import org.broad.igv.feature.xome.ExomeReferenceFrame;
import org.broad.igv.feature.xome.XomeUtils;
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
    static boolean exomeMode = false;

    static {
        //TODO This is a hack.
        if (!Globals.isHeadless()) {
            frames.add(getDefaultFrame());
        }
    }

    public static ReferenceFrame getDefaultFrame() {
        if (defaultFrame == null) {
            defaultFrame = new ReferenceFrame("genome");
        }
        return defaultFrame;
    }


    public static void setExomeMode(boolean b) {
        if(b == exomeMode) return;
        if(b) {
            switchToExomeMode();
        }
        else {
            switchToGenomeMode();
        }
    }


    public static boolean isExomeMode() {
        return exomeMode;
    }

    private static void switchToExomeMode() {
        List<Block> blocks = XomeUtils.getBlocks(defaultFrame.getChrName());
        ExomeReferenceFrame exomeFrame = new ExomeReferenceFrame(defaultFrame, blocks);
        defaultFrame = exomeFrame;
        frames.clear();
        frames.add(defaultFrame);
        exomeMode = true;
    }

    private static void switchToGenomeMode() {
        ReferenceFrame refFrame = new ReferenceFrame(defaultFrame);
        defaultFrame = refFrame;
        frames.clear();
        frames.add(defaultFrame);
        exomeMode = false;
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
                getDefaultFrame().setInterval(locus);
            }
        }
        frames.add(getDefaultFrame());
        getDefaultFrame().recordHistory();
    }


    public static void resetFrames(GeneList gl) {

        frames.clear();

        if (gl == null) {
            frames.add(getDefaultFrame());
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
                    getDefaultFrame().jumpTo(locus.getChr(), locus.getStart(), locus.getEnd());
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
                    Genome genome = GenomeManager.getInstance().getCurrentGenome();
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

