/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.util.MessageUtils;

/**
 * A class for performing search actions.  The class takes a view context and
 * search string as parameters.   The search string can be either
 * (a) a feature (e.g. gene),  or
 * (b) a locus string in the UCSC form,  e.g. chr1:100,000-200,000
 * <p/>
 * Note:  Currently the only recognized features are genes
 *
 * @author jrobinso
 */

/**
 * A class for performing search actions.  The class takes a view context and
 * search string as parameters.   The search string can be either
 * (a) a feature (e.g. gene),  or
 * (b) a locus string in the UCSC form,  e.g. chr1:100,000-200,000
 * <p/>
 * Note:  Currently the only recognized features are genes
 *
 * @author jrobinso
 */
public class SearchCommand implements Command {

    private static Logger log = Logger.getLogger(SearchCommand.class);

    String searchString;
    ReferenceFrame referenceFrame;
    boolean recordHistory = true;
    Genome genome;

    public SearchCommand(ReferenceFrame referenceFrame, String searchString) {
        this.referenceFrame = referenceFrame;
        this.searchString = searchString.trim();
        genome = GenomeManager.getInstance().getCurrentGenome();
    }

    public SearchCommand(ReferenceFrame referenceFrame, String searchString, boolean recordHistory) {
        this(referenceFrame, searchString);
        this.recordHistory = recordHistory;
    }


    public void execute() {

        if (log.isDebugEnabled()) {
            log.debug("Run search: " + searchString);
        }

        boolean success = false;

        // Space delimited?
        String[] tokens = searchString.split("\\s+");
        if (tokens.length >= 2) {
            String chr = tokens[0].trim();
            int start = Integer.parseInt(tokens[1].trim()) - 1; // Convert to UCSC convention
            int end = start + 1;
            if (tokens.length > 2) {
                end = Integer.parseInt(tokens[2].trim());
            }
            referenceFrame.jumpTo(chr, start, end);
            success = true;
        }

        // Feature search

        else {
            NamedFeature feature = FeatureDB.getFeature(searchString.toUpperCase().trim());
            if (feature != null) {
                int flankingRegion = PreferenceManager.getInstance().getAsInt(PreferenceManager.FLANKING_REGION);
                int start = Math.max(0, feature.getStart() - flankingRegion);
                int end = feature.getEnd() + flankingRegion;

                if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SEARCH_ZOOM)) {
                    referenceFrame.jumpTo(feature.getChr(), start, end);
                } else {
                    int center = (start + end) / 2;
                    referenceFrame.centerOnLocation(feature.getChr(), center);
                }

                if (log.isDebugEnabled()) {
                    log.debug("End search: " + searchString);
                }
                success = true;
            }


            // Apparently not a feature. Either a locus or track name.  Track names can be quoted,
            // loci are never quoted.
            else if (!searchString.contains("\"")) {
                String chr = null;
                int[] startEnd = null;
                int colon = searchString.indexOf(":");

                if (colon > 0) {

                    // The chromosome is that portion of the search string up to the colon.
                    chr = searchString.substring(0, colon);
                    String posString = searchString.substring(colon).replace(":", "");
                    startEnd = getStartEnd(posString);

                    if (startEnd != null) {
                        referenceFrame.jumpTo(chr, startEnd[0], startEnd[1]);

                        if (log.isDebugEnabled()) {
                            log.debug("End search: " + searchString);
                        }
                        success = true;
                    }
                } else {

                    // No chromosome delimiter (color),  The search string is either chromosome name
                    // or a locus in the current chromosome.
                    if (searchString.contains("-")) {

                        // Presense of a dash indicates this is a locus string in the current chromosome
                        startEnd = getStartEnd(searchString);
                        if (startEnd != null) {
                            referenceFrame.jumpTo(null, startEnd[0], startEnd[1]);

                            if (log.isDebugEnabled()) {
                                log.debug("End search: " + searchString);
                            }
                            success = true;

                        }
                    } else {

                        // No dash, this is either a chromosome or an unkown search string
                        Chromosome chromosome = GenomeManager.getInstance().getCurrentGenome().getChromosome(searchString);
                        if (chromosome != null || searchString.equals(Globals.CHR_ALL)) {
                            referenceFrame.setChromosomeName(searchString, true);
                            IGVMainFrame.getInstance().repaintDataAndHeaderPanels();
                            IGVMainFrame.getInstance().repaintStatusAndZoomSlider();

                            if (log.isDebugEnabled()) {
                                log.debug("End search: " + searchString);
                            }
                            success = true;
                        }

                    }
                }
            }
        }


        if (success) {
            if (recordHistory) {
                IGVMainFrame.getInstance().getSession().getHistory().push(searchString);
            }
        } else {
            if (!IGVMainFrame.getInstance().scrollToTrack(searchString.replaceAll("\"", ""))) {
                showError("Cannot find feature or locus: " + searchString);
            }
        }

        if (log.isDebugEnabled()) {
            log.debug("End search: " + searchString);
        }

    }

    /**
     * Return the start and end positions as a 2 element array for the input
     * position string.  UCSC conventions  are followed for coordinates,
     * specifically the internal representation is "zero" based (first base is
     * numbered 0) but the display representation is "one" based (first base is
     * numbered 1).   Consequently 1 is substracted from the parsed positions
     */
    private int[] getStartEnd(String posString) {
        try {
            String[] posTokens = posString.split("-");
            String startString = posTokens[0].replaceAll(",", "");
            int start = Math.max(0, Integer.parseInt(startString)) - 1;

            // Default value for end

            int end = start + 1;
            if (posTokens.length > 1) {
                String endString = posTokens[1].replaceAll(",", "");
                end = Integer.parseInt(endString);
            }

            if (posTokens.length == 1 || (end - start) < 10) {
                int center = (start + end) / 2;
                start = center - 20;
                end = center + 20;
            } else {
                String endString = posTokens[1].replaceAll(",", "");

                // Add 1 bp to end position t make it "inclusive"
                end = Integer.parseInt(endString);
            }

            return new int[]{Math.min(start, end), Math.max(start, end)};
        } catch (NumberFormatException numberFormatException) {
            return null;
        }

    }

    /**
     * TODO -- should this be here?  If not here where?
     *
     * @param message
     */
    private void showError(String message) {
        MessageUtils.showMessage(message);
    }
}
