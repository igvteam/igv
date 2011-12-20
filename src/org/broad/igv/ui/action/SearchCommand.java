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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.NamedFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.lists.GeneList;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;

import java.util.ArrayList;
import java.util.List;

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
        this(referenceFrame, searchString,
                IGV.getInstance().getGenomeManager().getCurrentGenome());
    }

    public SearchCommand(ReferenceFrame referenceFrame, String searchString, boolean recordHistory) {
        this(referenceFrame, searchString);
        this.recordHistory = recordHistory;
    }

    SearchCommand(ReferenceFrame referenceFrame, String searchString, Genome genome) {
        this.referenceFrame = referenceFrame;
        this.searchString = searchString.trim();
        this.genome = genome;
    }


    public void execute() {

        if (log.isDebugEnabled()) {
            log.debug("Run search: " + searchString);
        }

        SearchResult result = runSearch(searchString);

        if (result != null) {
            showSearchResult(result);
        } else {
            if (!IGV.getInstance().scrollToTrack(searchString.replaceAll("\"", ""))) {
                MessageUtils.showMessage("Cannot find feature or locus: " + searchString);
            }
        }

        if (log.isDebugEnabled()) {
            log.debug("End search: " + searchString);
        }
    }

    /**
     * Given a string, search for the appropriate data to show the user.
     * Differerent syntaxes are accepted.
     *
     * @param searchString
     * @return result
     *         SearchResult describing the results of the search. Will never
     *         be null, field type will equal SearchType.ERROR if something went wrong.
     */
    SearchResult runSearch(String searchString) {

        SearchResult result = new SearchResult(SearchType.ERROR, null, -1, -1, searchString);

        // Space delimited?
        String[] tokens = searchString.split("\\s+");
        if (tokens.length >= 2) {
            result = parseTokens(tokens);
            // Feature search
        } else {
            NamedFeature feature = FeatureDB.getFeature(searchString.toUpperCase().trim());
            if (feature != null) {
                return new SearchResult(feature, searchString);
            } else {
                //Map<String, NamedFeature> features = FeatureDB.getFeatures(searchString.toUpperCase().trim());
            }

            // Apparently not a feature. Either a locus or track name.  Track names can be quoted,
            // loci are never quoted.
            if (!searchString.contains("\"")) {
                result = calcChromoLocus(searchString);
            }
        }

        return result;
    }

    private void showSearchResult(SearchResult result) {
        int origZoom = referenceFrame.getZoom();
        if (result == null)
            result = new SearchResult(SearchType.ERROR, null, -1, -1, searchString);

        boolean showMessage = false;
        boolean success = true;
        String message = "Invalid search string: " + result.searchString;

        if (result.type != SearchType.ERROR && FrameManager.isGeneListMode()) {
            IGV.getInstance().getSession().setCurrentGeneList(null);
            IGV.getInstance().resetFrames();
        }

        switch (result.type) {
            case FEATURE:
                showFlankedRegion(result.chr, result.start, result.end);
                break;
            case LOCUS:
                referenceFrame.jumpTo(result.chr, result.start, result.end);
                break;
            case LOCI:
                message = result.getMessage();
                showMessage = message != null;
                GeneList geneList = new GeneList("", result.loci, false);
                IGV.getInstance().getSession().setCurrentGeneList(geneList);
                IGV.getInstance().resetFrames();
                break;
            case CHROMOSOME:
                referenceFrame.setChromosomeName(result.chr, true);
                IGV.getInstance().repaintDataAndHeaderPanels();
                IGV.getInstance().repaintStatusAndZoomSlider();
                break;
            case ERROR:
            default:
                success = false;
                showMessage = true;
        }
        if (success && recordHistory)
            IGV.getInstance().getSession().getHistory().push(searchString, origZoom);

        if (showMessage) {
            MessageUtils.showMessage(message);
        }

    }

    /**
     * Determine searchResult for white-space delimited search query.
     *
     * @param tokens
     * @return searchResult
     */
    private SearchResult parseTokens(String[] tokens) {
        SearchResult result;
        boolean success = false;
        int start = 0, end = 0;
        String chr = genome.getChromosomeAlias(tokens[0].trim());
        try {
            start = Integer.parseInt(tokens[1].trim()) - 1; // Convert to UCSC convention
            end = start + 1;
            if (tokens.length >= 3) {
                end = Integer.parseInt(tokens[2].trim());
            }
            result = new SearchResult(SearchType.LOCUS, chr, start, end, searchString);
            return result;
        } catch (NumberFormatException e) {
            // Multiple tokens, On the fly gene list ?
        }

        List<String> loci = new ArrayList<String>(tokens.length);
        Chromosome chromo = null;
        for (String t : tokens) {
            Locus l = new Locus(t);
            if (l.isValid()) {
                loci.add(t);
            } else if (FeatureDB.getFeature(t.toUpperCase().trim()) != null) {
                loci.add(t);
            } else if ((chromo = genome.getChromosome(t)) != null) {
                //Locus will only be valid if fully qualified.
                //Should rethink this approach
                String ft = t + ":1-" + chromo.getLength();
                Locus fl = new Locus(ft);
                if (fl.isValid())
                    loci.add(ft);
            }
        }
        result = new SearchResult(loci);

        if (loci.size() <= 1)
            return null;

        if (loci.size() != tokens.length) {
            result.setMessage("Not all portions of the search string are recognized as genes or loci.");
            //MessageUtils.showMessage("Not all portions of the search string are recognized as genes or loci.");
        }
        return result;
    }

    private SearchResult calcChromoLocus(String searchString) {
        int colonIdx = searchString.lastIndexOf(":");
        if (colonIdx > 0) {
            // The chromosome is that portion of the search string up to the last colon.
            String chr = genome.getChromosomeAlias(searchString.substring(0, colonIdx));
            String posString = searchString.substring(colonIdx).replace(":", "");
            int[] startEnd = getStartEnd(posString);
            if (startEnd != null) {
                return new SearchResult(SearchType.LOCUS, chr, startEnd[0], startEnd[1], searchString);
            }
        } else {
            // No chromosome delimiter (color),  The search string is either chromosome name
            // or a locus in the current chromosome.
            if (searchString.contains("-")) {
                // Presense of a dash indicates this is a locus string in the current chromosome
                int[] startEnd = getStartEnd(searchString);
                return new SearchResult(SearchType.LOCUS, null, startEnd[0], startEnd[1], searchString);
            } else {
                // No dash, this is either a chromosome or an unkown search string
                String chr = genome.getChromosomeAlias(searchString);
                Chromosome chromosome = genome.getChromosome(chr);
                if (chromosome != null || searchString.equals(Globals.CHR_ALL)) {
                    return new SearchResult(SearchType.CHROMOSOME, chr, 0, 0, searchString);
                }
            }
        }
        return new SearchResult(SearchType.ERROR, null, -1, -1, searchString);
    }

    private void showFlankedRegion(String chr, int start, int end) {
        int flankingRegion = PreferenceManager.getInstance().getAsInt(PreferenceManager.FLANKING_REGION);
        start = Math.max(0, start - flankingRegion);
        end = end + flankingRegion;

        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SEARCH_ZOOM)) {
            referenceFrame.jumpTo(chr, start, end);
        } else {
            int center = (start + end) / 2;
            referenceFrame.centerOnLocation(chr, center);
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

    enum SearchType {
        FEATURE,
        FEATURES,
        LOCUS,
        LOCI,
        TRACK,
        CHROMOSOME,
        ERROR
    }

    /*
    Container class for search results
     */
    class SearchResult {
        String chr;
        private int start;
        private int end;
        SearchType type;
        private String searchString;

        private List<String> loci;
        private String message;

        public SearchResult(SearchType type, String chr, int start, int end, String searchString) {
            this.type = type;
            this.chr = chr;
            this.start = start;
            this.end = end;
            this.searchString = searchString;
        }

        public SearchResult(List<String> loci) {
            this.type = SearchType.LOCI;
            this.loci = loci;
        }

        public SearchResult(NamedFeature feature, String searchString) {
            this(SearchType.FEATURE, feature.getChr(), feature.getStart(), feature.getEnd(),
                    searchString);
        }

        private void setMessage(String message) {
            this.message = message;
        }

        public String getMessage() {
            return this.message;
        }

        List<String> getLoci() {
            return this.loci;
        }
    }
}