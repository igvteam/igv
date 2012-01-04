/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTIES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */


package org.broad.igv.ui.action;

//~--- non-JDK imports --------------------------------------------------------

import com.sun.org.apache.regexp.internal.RE;
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

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
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
    public static int SEARCH_LIMIT = 20;
    private boolean askUser = false;

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

        List<SearchResult> results = runSearch(searchString);
        if (askUser) {
            results = askUserFeature(results);
            if (results == null) {
                if (log.isDebugEnabled()) {
                    log.debug("Multiple results, show cancelled: " + searchString);
                }
                return;
            }
        }

        showSearchResult(results);

        if (log.isDebugEnabled()) {
            log.debug("End search: " + searchString);
        }
    }

    /**
     * Given a string, search for the appropriate data to show the user.
     * Different syntaxes are accepted.
     * <p/>
     * In general, whitespace delimited tokens are treated separately and each are shown.
     * There is 1 exception to this. A locus of form chr1   1   10000 will be treated the same
     * as chr1:1-10000. Only one entry of this form can be entered, chr1    1   10000 chr2:1-1000 will
     * not be recognized.
     *
     * @param searchString Feature name (EGFR), chromosome (chr1), or locus string (chr1:1-100 or chr1:6)
     *                     Partial matches to a feature name (EG) will return multiple results, and
     *                     ask the user which they want.
     * @return result
     *         List<SearchResult> describing the results of the search. Will never
     *         be null, field type will equal ResultType.ERROR if something went wrong.
     */
    List<SearchResult> runSearch(String searchString) {

        List<SearchResult> results = new ArrayList<SearchResult>();

        searchString = searchString.replace("\"", "");

        ResultType wholeStringType = checkTokenType(searchString);
        if (wholeStringType == ResultType.LOCUS) {
            results.add(calcChromoLocus(searchString));
            return results;
        }

        // Space delimited?
        String[] tokens = searchString.split("\\s+");
        for (String s : tokens) {
            results.addAll(parseToken(s));
        }

        if (results.size() == 0) {
            SearchResult result = new SearchResult();
            result.setMessage("Invalid Search String: " + searchString);
            results.add(result);
        }

        return results;
    }

    public void showSearchResult(List<SearchResult> results) {
        int origZoom = referenceFrame.getZoom();
        SearchResult result = new SearchResult();
        if (results == null || results.size() == 0) {
            results = new ArrayList<SearchResult>();
            results.add(result);
        }
        boolean showMessage = false;
        boolean success = true;
        String message = "Invalid search string: " + searchString;

        if (results.size() == 1) {
            result = results.get(0);
            if (result.type != ResultType.ERROR && FrameManager.isGeneListMode()) {
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
                case CHROMOSOME:
                    referenceFrame.setChromosomeName(result.chr, true);
                    IGV.getInstance().repaintDataAndHeaderPanels();
                    IGV.getInstance().repaintStatusAndZoomSlider();
                    break;
                case ERROR:
                default: {
                    message = "Cannot find feature or locus: " + searchString;
                    success = false;
                    showMessage = true;
                }

            }
        } else {
            List<String> loci = new ArrayList<String>(results.size());
            message = "";
            for (SearchResult res : results) {
                if (res.type != ResultType.ERROR) {
                    loci.add(res.getLocus());
                } else {
                    message = message + res.getMessage() + "\n";
                    showMessage = true;
                }
            }
            GeneList geneList = new GeneList("", loci, false);
            IGV.getInstance().getSession().setCurrentGeneList(geneList);
            IGV.getInstance().resetFrames();

        }


        if (success && recordHistory) {
            IGV.getInstance().getSession().getHistory().push(searchString, origZoom);
        }
        if (showMessage) {
            MessageUtils.showMessage(message);
        }

    }

    /**
     * Get a list of strings of feature names suitable for display, containing only
     * those search results which were not an error
     *
     * @param results
     * @param longName Whether to use the long (true) or short (false)
     *                 of search results.
     * @return Array of strings of results found.
     */
    public static Object[] getSelectionList(List<SearchResult> results, boolean longName) {
        ArrayList<String> options = new ArrayList<String>(Math.min(results.size(), SEARCH_LIMIT));
        for (SearchResult result : results) {
            if (result.type == ResultType.ERROR) {
                continue;
            }
            if (longName) {
                options.add(result.getLongName());
            } else
                options.add(result.getShortName());
        }

        return options.toArray();
    }

    /**
     * Display a dialog asking user which search result they want
     * to display. Number of results are limited to SEARCH_LIMIT.
     * The user can select multiple options, in which case all
     * are displayed.
     *
     * @param results
     * @return SearchResults which the user has selected.
     *         Will be null if cancelled
     */
    private List<SearchResult> askUserFeature(List<SearchResult> results) {

        Object[] list = getSelectionList(results, true);
        JList ls = new JList(list);
        ls.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        //ls.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);

        final JOptionPane pane = new JOptionPane(ls, JOptionPane.PLAIN_MESSAGE, JOptionPane.OK_CANCEL_OPTION);
        final Dialog dialog = pane.createDialog("Features");
        dialog.setModalityType(Dialog.ModalityType.APPLICATION_MODAL);

        //On double click, show that option
        ls.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() >= 2) {
                    dialog.setVisible(false);
                    pane.setValue(JOptionPane.OK_OPTION);
                    dialog.dispose();
                }
            }
        });

        dialog.setVisible(true);

        int resp = (Integer) pane.getValue();

        List<SearchResult> val = null;
        if (resp == JOptionPane.OK_OPTION) {
            int[] selected = ls.getSelectedIndices();
            val = new ArrayList<SearchResult>(selected.length);
            for (int ii = 0; ii < selected.length; ii++) {
                val.add(ii, results.get(selected[ii]));
            }
        }
        return val;

    }

    /**
     * Check token type using regex.
     *
     * @param token
     * @return
     */
    ResultType checkTokenType(String token) {
        //Regexp for a number with commas in it (no periods)
        String num_withcommas = "(((\\d)+,?)+)";
        //Not ideal, will match chr + (1 or 2 digit number) or chr[X,Y, or M]
        String chromo_string = "^chr([\\d]{1,2}|[XYM])";
        RE chromo = new RE(chromo_string + "\\s*$", RE.MATCH_CASEINDEPENDENT);
        //This will match chr1:1-100, chr1:1, chr1  1, chr1 1   100
        RE chromo_range = new RE(chromo_string + "(:|(\\s)+)" + num_withcommas + "(-|(\\s)+)?" + num_withcommas + "?(\\s)*$",
                RE.MATCH_CASEINDEPENDENT);

        //Simple feature, which is letters/numbers only
        RE feature = new RE("^(\\s)*(\\w)+(\\s)*$", RE.MATCH_CASEINDEPENDENT);
        //Mutation notation. e.g. KRAS:G12C
        RE feature_mut = new RE("^(\\w)+:[A-Z]" + num_withcommas + "[A-Z](\\s)*$", RE.MATCH_CASEINDEPENDENT);
        if (chromo.match(token)) {
            return ResultType.CHROMOSOME;
        } else if (chromo_range.match(token)) {
            return ResultType.LOCUS;
        } else if (feature.match(token)) {
            return ResultType.FEATURE;
        } else if (feature_mut.match(token)) {
            return ResultType.FEATURE_MUT;
        }

        return ResultType.ERROR;
    }

    /**
     * Determine searchResult for white-space delimited search query.
     *
     * @param token
     * @return searchResult
     */
    private List<SearchResult> parseToken(String token) {

        List<SearchResult> results = new ArrayList<SearchResult>();
        List<NamedFeature> features;

        //Guess at token type via regex.
        //We don't assume success
        ResultType type = checkTokenType(token);
        SearchResult result;
        switch (type) {
            case LOCUS:
                //Check if a full or partial locus string
                result = calcChromoLocus(token);
                if (result.type != ResultType.ERROR) {
                    results.add(result);
                    return results;
                }
                break;
            case FEATURE_MUT:
                //We know it has the right form, but may
                //not be valid feature name or mutation
                //which exists.
                String[] items = token.split(":");
                String name = items[0].trim().toUpperCase();
                String coords = items[1];
                char refAASymbol = coords.subSequence(0, 1).charAt(0);
                coords = coords.substring(1, coords.length() - 1);
                int location = Integer.parseInt(coords) - 1;
                features = FeatureDB.getMutation(name, location, refAASymbol);
                //Only keep the largest one
                if (features.size() >= 1) {
                    NamedFeature feat = features.get(0);
                    //result = new SearchResult(feat);
                    //Zoom in on mutation of interest
                    result = new SearchResult(ResultType.LOCUS, feat.getChr(), Math.max(0, location));
                    results.add(result);
                    return results;
                }
                break;
            case CHROMOSOME:
                String chr = genome.getChromosomeAlias(token);
                Chromosome chromo = genome.getChromosome(chr);
                if (chromo != null) {
                    result = new SearchResult(ResultType.CHROMOSOME, chr, 1, chromo.getLength());
                    results.add(result);
                    return results;
                }
                //Chromosome string can look similar to feature.
                //If we fail at chromosome, check against being a feature
            default:
            case FEATURE:
                //Check if a feature
                NamedFeature feat = FeatureDB.getFeature(token.toUpperCase().trim());
                if (feat != null) {
                    results.add(new SearchResult(feat));
                    return results;
                }

                //Check inexact match
                //We will later want to ask the user which of these to keep
                features = FeatureDB.getFeaturesList(searchString, SEARCH_LIMIT);
                if (features.size() > 0) {
                    askUser |= features.size() >= 2;
                    return getResults(features);
                }
        }

        result = new SearchResult();
        result.setMessage("Invalid token: " + token);
        results.add(result);
        return results;

    }

    /**
     * Parse a string of locus coordinates.
     * Can have whitespace delimiters, and be mising second coordinate,
     * but must have 1st coordinate.
     *
     * @param searchString
     * @return
     */
    private SearchResult calcChromoLocus(String searchString) {
        int colonIdx = searchString.lastIndexOf(":");
        int[] startEnd = null;
        String chr = null;
        if (colonIdx > 0) {
            chr = searchString.substring(0, colonIdx);
            String posString = searchString.substring(colonIdx).replace(":", "");
            startEnd = getStartEnd(posString);
        } else {
            //Assume whitespace delimited
            String[] tokens = searchString.split("\\s+");
            chr = tokens[0];
            String posString = tokens[1];
            if (tokens.length >= 3) {
                posString += "-" + tokens[2];
            }
            startEnd = getStartEnd(posString);

        }
        chr = genome.getChromosomeAlias(chr);
        Chromosome chromosome = genome.getChromosome(chr);
        if (chromosome != null && !searchString.equals(Globals.CHR_ALL)) {
            if (startEnd != null) {
                return new SearchResult(ResultType.LOCUS, chr, startEnd[0], startEnd[1]);
            }
            return new SearchResult(ResultType.CHROMOSOME, chr, 0, 0);
        }
        return new SearchResult(ResultType.ERROR, chr, -1, -1);
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

    enum ResultType {
        FEATURE,
        FEATURE_MUT,
        LOCUS,
        CHROMOSOME,
        ERROR
    }

    /*
    Container class for search results
     */
    public static class SearchResult {
        String chr;
        private int start;
        private int end;
        ResultType type;

        private String locus;
        private String message;
        private NamedFeature feature;
        private String coords;

        public SearchResult() {
            this(ResultType.ERROR, null, -1, -1);
        }

        public SearchResult(ResultType type, String chr, int start, int end) {
            this.type = type;
            this.chr = chr;
            this.start = start;
            this.end = end;
            this.coords = Locus.getFormattedLocusString(chr, start + 1, end);
            this.locus = this.coords;
        }

        public SearchResult(NamedFeature feature) {
            this(ResultType.FEATURE, feature.getChr(), feature.getStart(), feature.getEnd());
            this.feature = feature;
            this.locus = this.feature.getName();
        }

        public SearchResult(ResultType type, String chr, int location) {
            this.type = type;
            this.chr = chr;
            this.start = location;
            this.end = location + 1;
            this.locus = chr + ":" + (location + 1);
        }

        void setMessage(String message) {
            this.message = message;
        }

        public String getMessage() {
            return this.message;
        }

        /**
         * Always a coordinate string.
         * eg chr1:1-100
         *
         * @return
         */
        private String getCoordinates() {
            return this.coords;
        }

        /**
         * Either a feature name, or coordinates
         *
         * @return
         */
        String getLocus() {
            return this.locus;
        }

        String getShortName() {
            if (this.type == ResultType.FEATURE) {
                return this.feature.getName();
            } else {
                return this.getLocus();
            }
        }

        /**
         * Format for display. If a feature,
         * Featurename (chromosome:start-end)
         * eg EGFR (chr7:55,054,218-55,242,525)
         * <p/>
         * Otherwise, just locus
         *
         * @return
         */
        String getLongName() {
            if (this.type == ResultType.FEATURE) {
                return feature.getName() + " (" + this.getCoordinates() + ")";
            } else {
                return this.getLocus();
            }
        }

    }

    /**
     * Get a list of search results from the provided objects,
     * which must be loci strings or NamedFeature objects.
     *
     * @param objects
     * @return
     */
    public static List<SearchResult> getResults(List<NamedFeature> objects) {
        List<SearchResult> results = new ArrayList<SearchResult>(objects.size());
        for (NamedFeature f : objects) {
            results.add(new SearchCommand.SearchResult(f));
        }
        return results;
    }


}