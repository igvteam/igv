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


package org.broad.igv.ui.action;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.lists.GeneList;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.tribble.Feature;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.*;
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
        this(referenceFrame, searchString, GenomeManager.getInstance().getCurrentGenome());
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
    public List<SearchResult> runSearch(String searchString) {

        List<SearchResult> results = new ArrayList<SearchResult>();

        searchString = searchString.replace("\"", "");

        Set<ResultType> wholeStringType = checkTokenType(searchString);
        if (wholeStringType.contains(ResultType.LOCUS)) {
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
            if (result.type != ResultType.ERROR) {//FrameManager.isGeneListMode()) {
                IGV.getInstance().getSession().setCurrentGeneList(null);
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
                    IGV.repaintPanelsHeadlessSafe();
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
        }

        IGV.getInstance().resetFrames();


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
     * Intended to be inclusive, returns all possible matches
     *
     * @param token
     * @return
     */
    Set<ResultType> checkTokenType(String token) {
        token = token.trim();

        //Regexp for a number with commas in it (no periods)
        String num_withcommas = "(((\\d)+,?)+)";

        //chromosome can include anything except whitespace
        String chromo_string = "(\\S)+";


        String chromo = chromo_string;
        //This will match chr1:1-100, chr1:1, chr1  1, chr1 1   100
        String chromo_range = chromo_string + "(:|(\\s)+)" + num_withcommas + "(-|(\\s)+)?" + num_withcommas + "?(\\s)*";

        //Simple feature
        String feature = chromo_string;
        //Amino acid mutation notation. e.g. KRAS:G12C. * is stop codon
        String featureMutAA = chromo_string + ":[A-Z,a-z,*]" + num_withcommas + "[A-Z,a-z,*]";

        //Nucleotide mutation notation. e.g. KRAS:123A>T
        String nts = "[A,C,G,T,a,c,g,t]";
        String featureMutNT = chromo_string + ":" + num_withcommas + nts + "\\>" + nts;

        Set<ResultType> possibles = new HashSet<ResultType>();
        Map<ResultType, String> matchers = new HashMap<ResultType, String>();
        matchers.put(ResultType.CHROMOSOME, chromo);
        matchers.put(ResultType.FEATURE, feature);
        matchers.put(ResultType.LOCUS, chromo_range);
        matchers.put(ResultType.FEATURE_MUT_AA, featureMutAA);
        matchers.put(ResultType.FEATURE_MUT_NT, featureMutNT);
        for (ResultType type : matchers.keySet()) {
            if (token.matches(matchers.get(type))) { //note: entire string must match
                possibles.add(type);
            }
        }

        return possibles;
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
        Set<ResultType> types = checkTokenType(token);
        SearchResult result;
        if (types.contains(ResultType.LOCUS) || types.contains(ResultType.CHROMOSOME)) {
            //Check if a full or partial locus string
            result = calcChromoLocus(token);
            if (result.type != ResultType.ERROR) {
                results.add(result);
                return results;
            }
        }

        //2 possible mutation notations, either amino acid (A123B) or nucleotide (123G>C)
        if (types.contains(ResultType.FEATURE_MUT_AA) || types.contains(ResultType.FEATURE_MUT_NT)) {
            //We know it has the right form, but may
            //not be valid feature name or mutation
            //which exists.
            String[] items = token.toUpperCase().split(":");
            String name = items[0].trim().toUpperCase();
            String coords = items[1];
            int coordLength = coords.length();

            Map<Integer, BasicFeature> genomePosList;

            //Should never match both mutation notations
            if (types.contains(ResultType.FEATURE_MUT_AA)) {
                String refSymbol = coords.substring(0, 1);
                String mutSymbol = coords.substring(coordLength - 1);

                String strLoc = coords.substring(1, coordLength - 1);
                int location = Integer.parseInt(strLoc) - 1;

                genomePosList = FeatureDB.getMutationAA(name, location + 1, refSymbol, mutSymbol, genome);
            } else if (types.contains(ResultType.FEATURE_MUT_NT)) {
                //Exclude the "A>T" at end
                String strLoc = coords.substring(0, coordLength - 3);
                String refSymbol = coords.substring(coordLength - 3, coordLength - 2);
                int location = Integer.parseInt(strLoc) - 1;
                genomePosList = FeatureDB.getMutationNT(name, location + 1, refSymbol, genome);
            } else {
                //This should never happen
                throw new IllegalArgumentException("Something went wrong parsing input token");
            }
            askUser |= genomePosList.size() >= 2;

            for (int genomePos : genomePosList.keySet()) {
                Feature feat = genomePosList.get(genomePos);
                //Zoom in on mutation of interest
                //The +2 accounts for centering on the center of the amino acid, not beginning
                //and converting from 0-based to 1-based (which getStartEnd expects)
                int[] locs = getStartEnd("" + (genomePos + 2));
                result = new SearchResult(ResultType.LOCUS, feat.getChr(), locs[0], locs[1]);
                results.add(result);

            }
            return results;
        }

        if (types.contains(ResultType.FEATURE)) {
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
     * Can have whitespace delimiters, and be missing second coordinate,
     * but must have 1st coordinate.
     *
     * @param searchString
     * @return
     */
    private SearchResult calcChromoLocus(String searchString) {
        /*
        chromosome can have whitespace or : delimiter
        chromosome also might have : in the name
         */
        int[] startEnd = null;
        String[] tokens = searchString.split("\\s+");

        String chr = tokens[0];
        boolean whitespace_delim = tokens.length >= 2;
        if (whitespace_delim) {
            String posString = tokens[1];
            if (tokens.length >= 3) {
                posString += "-" + tokens[2];
            }
            startEnd = getStartEnd(posString);
        } else {
            //Not whitespace delimited
            //Could be chromoname:1-100, chromoname:1, chromoname

            int colonIdx = searchString.lastIndexOf(":");
            if (colonIdx > 0) {
                chr = searchString.substring(0, colonIdx);
                String posString = searchString.substring(colonIdx).replace(":", "");
                startEnd = getStartEnd(posString);
                //This MAY for case of chromoname having semicolon in it
                if (startEnd == null) {
                    chr = searchString;
                }
            }
        }

        //startEnd will have coordinates if found.
        chr = genome.getChromosomeAlias(chr);
        Chromosome chromosome = genome.getChromosome(chr);
        //If we couldn't find chromosome, check
        //whole string
        if (chromosome == null) {
            chr = genome.getChromosomeAlias(tokens[0]);
            chromosome = genome.getChromosome(chr);
            if (chromosome != null) {
                //Found chromosome
                startEnd = null;
            }
        }

        if (chromosome != null && !searchString.equals(Globals.CHR_ALL)) {
            if (startEnd != null) {
                return new SearchResult(ResultType.LOCUS, chr, startEnd[0], startEnd[1]);
            }
            return new SearchResult(ResultType.CHROMOSOME, chr, 0, chromosome.getLength() - 1);
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
     * numbered 0) and end-exclusive, but the display representation is "one" based (first base is
     * numbered 1) and end-inclusive.   Consequently 1 is subtracted from the parsed positions
     */
    private static int[] getStartEnd(String posString) {
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
                int widen = 20;
                start = center - widen;
                start = Math.max(0, start);
                end = center + widen + 1;
            }

            return new int[]{Math.min(start, end), Math.max(start, end)};
        } catch (NumberFormatException numberFormatException) {
            return null;
        }

    }

    public enum ResultType {
        FEATURE,
        FEATURE_MUT_AA,
        FEATURE_MUT_NT,
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

        public ResultType getType() {
            return type;
        }

        public String getChr() {
            return chr;
        }

        public int getStart() {
            return start;
        }

        public int getEnd() {
            return end;
        }
    }

    /**
     * Get a list of search results from the provided objects,
     * which must be NamedFeature objects.
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