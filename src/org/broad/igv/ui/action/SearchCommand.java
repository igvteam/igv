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
    private static int SEARCH_LIMIT = 20;
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
        }

        showSearchResult(results);

        if (log.isDebugEnabled()) {
            log.debug("End search: " + searchString);
        }
    }

    /**
     * Given a string, search for the appropriate data to show the user.
     * Different syntaxes are accepted.
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
        // Space delimited?
        String[] tokens = searchString.split("\\s+");
        if (tokens.length >= 2) {
            results.addAll(parseTokens(tokens));
            // Feature search
        } else {
            //Check exact match
            NamedFeature feature = FeatureDB.getFeature(searchString.toUpperCase().trim());
            if (feature != null) {
                results.add(new SearchResult(feature));
                return results;
            } else {
                //Check inexact match
                //We will later want to ask the user which of these to keep
                List<NamedFeature> features = FeatureDB.getFeaturesList(searchString, SEARCH_LIMIT);
                if (features.size() > 0) {
                    askUser = features.size() >= 2;
                    return getResults(features);
                }
            }

            // Apparently not a feature.
            // Either a locus or track name.  Track names can be quoted,
            // loci are never quoted.

            //Update: not supporting track names here
            results.add(calcChromoLocus(searchString));
        }

        if (results.size() == 0) {
            SearchResult result = new SearchResult();
            result.setMessage("Invalid Search String: " + searchString);
            results.add(result);
        }

        return results;
    }

    private void showSearchResult(List<SearchResult> results) {
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
     * Display a dialog asking user which search result they want
     * to display. Number of results are limited to SEARCH_LIMIT.
     * The user can select multiple options, in which case all
     * are displayed.
     *
     * @param results
     * @return SearchResults which the user has selected
     */
    private List<SearchResult> askUserFeature(List<SearchResult> results) {

        String[] options = new String[Math.min(results.size(), SEARCH_LIMIT)];
        for (int ii = 0; ii < options.length; ii++) {
            options[ii] = results.get(ii).getLongName();
        }

        JList ls = new JList(options);
        ls.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

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
     * Tokens will be separated somehow. Differentiating token types may be easier
     * using regex.
     * @param token
     * @return
     */
/*    private TokenType parseIndividualToken(String token){
        //Regexp for a number with commas in it (no periods)
        String num_withcommas = "((\\d)+,?)+";
        String chromo_string = "^chr([\\d]{1,2}|[XY])";
        RE chromo = new RE(chromo_string + "\\s*$", RE.MATCH_CASEINDEPENDENT);
        RE chromo_range = new RE(chromo_string + ":|\\s+" + num_withcommas + "(-|\\s+)" + num_withcommas + "(\\s)*$");
        if(chromo.match(token)){
            return TokenType.CHROMOSOME;
        }else if (chromo_range.match(token)){
            return TokenType.LOCUS;
        }
        
        return TokenType.ERROR;
    }*/

    /**
     * Determine searchResult for white-space delimited search query.
     *
     * @param tokens
     * @return searchResult
     */
    private List<SearchResult> parseTokens(String[] tokens) {
        List<SearchResult> results = new ArrayList<SearchResult>();
        int start, end;
        String chr = genome.getChromosomeAlias(tokens[0].trim());
        try {
            start = Integer.parseInt(tokens[1].trim()) - 1; // Convert to UCSC convention
            end = start + 1;
            if (tokens.length >= 3) {
                end = Integer.parseInt(tokens[2].trim());
            }
            SearchResult result = new SearchResult(ResultType.LOCUS, chr, start, end);
            results.add(result);
            return results;
        } catch (NumberFormatException e) {
            // Multiple tokens, On the fly gene list ?
        }

        //List<String> loci = new ArrayList<String>(tokens.length);
        results = new ArrayList<SearchResult>(tokens.length);
        SearchResult result;
        Chromosome chromo;
        for (String t : tokens) {
            Locus l = new Locus(t);
            result = new SearchResult();
            result.setMessage("Invalid token: " + t);
            if (l.isValid()) {
                result = new SearchResult(l);
            } else if ((chromo = genome.getChromosome(t)) != null) {
                //Locus will only be valid if fully qualified.
                //Should rethink this approach
                String ft = t + ":1-" + chromo.getLength();
                Locus fl = new Locus(ft);
                if (fl.isValid()) {
                    result = new SearchResult(fl);
                }
            } else {
                NamedFeature feat = FeatureDB.getFeature(t.toUpperCase().trim());
                if (feat != null) {
                    result = new SearchResult(feat);
                }
            }
            results.add(result);
        }

        return results;
    }

    private SearchResult calcChromoLocus(String searchString) {
        int colonIdx = searchString.lastIndexOf(":");
        if (colonIdx > 0) {
            // The chromosome is that portion of the search string up to the last colon.
            String chr = genome.getChromosomeAlias(searchString.substring(0, colonIdx));
            String posString = searchString.substring(colonIdx).replace(":", "");
            int[] startEnd = getStartEnd(posString);
            if (startEnd != null) {
                return new SearchResult(ResultType.LOCUS, chr, startEnd[0], startEnd[1]);
            }
        } else {
            // No chromosome delimiter (color),  The search string is either chromosome name
            // or a locus in the current chromosome.
            if (searchString.contains("-")) {
                // Presence of a dash indicates this is a locus string in the current chromosome
                int[] startEnd = getStartEnd(searchString);
                return new SearchResult(ResultType.LOCUS, null, startEnd[0], startEnd[1]);
            } else {
                // No dash, this is either a chromosome or an unkown search string
                String chr = genome.getChromosomeAlias(searchString);
                Chromosome chromosome = genome.getChromosome(chr);
                if (chromosome != null || searchString.equals(Globals.CHR_ALL)) {
                    return new SearchResult(ResultType.CHROMOSOME, chr, 0, 0);
                }
            }
        }
        return new SearchResult(ResultType.ERROR, null, -1, -1);
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
        LOCUS,
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
            this.coords = Locus.getFormattedLocusString(chr, start, end);
            this.locus = this.coords;
        }

        public SearchResult(Locus locus) {
            this(ResultType.LOCUS, locus.getChr(), locus.getStart(), locus.getEnd());
        }

        public SearchResult(NamedFeature feature) {
            this(ResultType.FEATURE, feature.getChr(), feature.getStart(), feature.getEnd());
            this.feature = feature;
            this.locus = this.feature.getName();
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
         * <Featurename> (<chromosome>:<start>-<end>)
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
    List<SearchResult> getResults(List<NamedFeature> objects) {
        List<SearchResult> results = new ArrayList<SearchResult>(objects.size());
        for (NamedFeature f : objects) {
            results.add(new SearchResult(f));
        }
        return results;
    }


}