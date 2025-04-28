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


package org.broad.igv.ui.action;

//~--- non-JDK imports --------------------------------------------------------

import htsjdk.tribble.Feature;
import htsjdk.tribble.NamedFeature;
import org.broad.igv.logging.*;
import org.broad.igv.Globals;
import org.broad.igv.annotations.ForTesting;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.lists.GeneList;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.Track;
import org.broad.igv.ucsc.SearchAPI;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;

import java.net.URL;
import java.util.List;
import java.util.*;


/**
 * A class for performing search actions.  The class takes a view context and
 * search string as parameters.   The search string can be either
 * (a) a feature (e.g. gene),  or
 * (b) a locus string in the UCSC form,  e.g. chr1:100,000-200,000
 *
 * @author jrobinso
 */
public class SearchCommand implements Runnable {

    private static Logger log = LogManager.getLogger(SearchCommand.class);
    public static int SEARCH_LIMIT = 10000;

    String searchString;
    ReferenceFrame referenceFrame;
    boolean recordHistory = true;
    Genome genome;


    private static HashMap<String, ResultType> tokenMatchers;

    static String featureMutAA = "(\\S)+" + ":[A-Z,a-z,*]" + "(((\\d)+,?)+)" + "[A-Z,a-z,*]";
    static String featureMutNT = "(\\S)+" + ":" + "(\\S)+" + "[A,C,G,T,a,c,g,t]" + "\\>" + "[A,C,G,T,a,c,g,t]";


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


    public void run() {
        List<SearchResult> results = runSearch(searchString);
        showSearchResult(results);
    }

    /**
     * Given a string, search for the appropriate locus or loci.
     * Different syntaxes are accepted.
     * <p/>
     * In general, space delimited tokens are treated separately and each are shown.
     * There is 1 exception to this. A locus of form chr1   1   10000 will be treated the same
     * as chr1:1-10000. Only one entry of this form can be entered, chr1    1   10000 chr2:1-1000 will
     * not be recognized.
     *
     * @param searchString Feature name (EGFR), chromosome (chr1), or locus string (chr1:1-100 or chr1:6)
     *                     Partial matches to a feature name (EG) will return multiple results, and
     *                     ask the user which they want.
     * @return result
     * List<SearchResult> describing the results of the search. Will never
     * be null, field type will equal ResultType.ERROR if something went wrong.
     */
    public List<SearchResult> runSearch(String searchString) {

        // Check for special "liftover" syntax.  This allows searching based on coordinates from another genome
        // (the "target" genome) if an associated liftover map is defined for the target genome.
//        Liftover liftover = null;
//        if (searchString.startsWith("!") && genome.getLiftoverMap() != null) {
//            int idx = searchString.indexOf(' ');
//            String genomeKey = searchString.substring(1, idx);
//            liftover = genome.getLiftoverMap().get(genomeKey);
//            if (liftover != null) {
//                searchString = searchString.substring(idx + 1);
//            }
//        }

        List<SearchResult> results = new ArrayList<>();

        searchString = searchString.replace("\"", "");

        // If the search string is space delimited see if it looks like a space delimited locus string (e.g. chr 100 200)
        String[] tokens = searchString.split("\\s+");
        if (tokens.length > 1 && tokens.length <= 3) {
            boolean mightBeLocus = true;
            for (int i = 1; i < tokens.length; i++) {
                mightBeLocus = mightBeLocus && isInteger(tokens[i]);
            }
            if (mightBeLocus) {
                Chromosome c1 = genome.getChromosome(tokens[0]);
                if (c1 != null) {
                    Chromosome c2 = genome.getChromosome(tokens[1]);
                    if (c2 == null) {
                        results.add(calcChromoLocus(searchString));
                        return results;
                    }
                }
            }
        }

        for (String s : tokens) {
            SearchResult result = parseToken(s);
            if (result != null) {
                results.add(result);
            }
        }


        // If this is a liftover search map the results
//        if (liftover != null) {
//            List<SearchResult> mappedResults = new ArrayList<>();
//            for (SearchResult result : results) {
//                if (result.getType() == ResultType.LOCUS) {
//                    List<Range> mapped = liftover.map(new Range(result.getChr(), result.getStart(), result.getEnd()));
//                    for (Range m : mapped) {
//                        mappedResults.add(new SearchResult(result.type, m.chr, m.start, m.end));
//                    }
//                } else {
//                    // ??? Error
//                }
//            }
//            results = mappedResults;
//        }

        return results;
    }

    public void showSearchResult(List<SearchResult> results) {
        int origZoom = referenceFrame.getZoom();
        if (results == null || results.size() == 0) {
            results = new ArrayList<>();
            results.add(new SearchResult());
        }

        boolean showMessage = false;
        boolean success = true;
        String message = "Invalid search string: " + searchString;

        boolean isGeneListMode = FrameManager.isGeneListMode();
        boolean resetFrames = false;

        if (results.size() == 1) {
            resetFrames = isGeneListMode;  // From gene list -> single locus

            SearchResult result = results.get(0);
            if (result.type != ResultType.ERROR) {//FrameManager.isGeneListMode()) {
                IGV.getInstance().getSession().setCurrentGeneList(null);
            }

            switch (result.type) {
                case FEATURE:
                    showFlankedRegion(result.chr, result.start, result.end);
                    break;
                case LOCUS:
                    if (result.chr.equalsIgnoreCase(Globals.CHR_ALL)) {
                        referenceFrame.changeChromosome(Globals.CHR_ALL, false);
                    } else {
                        Chromosome chromosome = GenomeManager.getInstance().getCurrentGenome().getChromosome(result.chr);
                        if (chromosome == null) {
                            message = "Unknow chromosome: " + result.chr;
                            success = false;
                            showMessage = true;
                        } else if (result.start > chromosome.getLength()) {
                            message = "Range " + result.locus + " is beyond the end of the chromosome";
                            success = false;
                            showMessage = true;
                        } else {
                            referenceFrame.jumpTo(result.chr, result.start, result.end);
                        }
                    }
                    break;
                case ERROR:
                default: {
                    message = "Cannot find feature or locus: " + searchString;
                    success = false;
                    showMessage = true;
                }
            }
        } else {
            resetFrames = true;  // New set of loci
            List<String> loci = new ArrayList<>(results.size());
            message = "<html>";
            for (SearchResult res : results) {
                if (res.type != ResultType.ERROR) {
                    loci.add(res.getLocus());
                } else {
                    message = message + res.getMessage() + "<br>";
                    showMessage = true;
                }
            }
            GeneList geneList = new GeneList("", loci);
            IGV.getInstance().getSession().setCurrentGeneList(geneList);
        }
        if (resetFrames) {
            IGV.getInstance().resetFrames();
        } else {
            IGV.getInstance().repaint();
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
     * Check token type using regex.
     * Intended to be inclusive, returns all possible matches
     *
     * @param token
     * @return
     */
    Set<ResultType> checkTokenType(String token) {
        token = token.trim();
        Set<ResultType> possibles = new HashSet<>();
        for (String key : tokenMatchers.keySet()) {
            if (token.matches(key)) {
                possibles.add(tokenMatchers.get(key));
            }
        }
        return possibles;
    }

    /**
     * Determine searchResult for string token.
     *
     * @param token
     * @return searchResult
     */
    private SearchResult parseToken(String token) {

        // Check featureDB first -- this is cheap
        NamedFeature feat = FeatureDB.getFeature(token.toUpperCase().trim());
        if (feat != null) {
            return new SearchResult(feat);
        }

        //Check if a full or partial locus string
        SearchResult result = calcChromoLocus(token);
        if (result != null) {
            return result;
        }

        //Check if we have an exact match for the feature name
        List<Track> searchableTracks = IGV.getInstance().getAllTracks().stream().filter(Track::isSearchable).toList();
        for (Track t : searchableTracks) {
            NamedFeature match = t.search(token);
            if (match != null) {
                return new SearchResult(match);
            }
        }

        // Try the webservice
        feat = searchWebservice(token);
        if (feat != null) {
            return new SearchResult(feat);
        }


        //2 possible mutation notations, either amino acid (A123B) or nucleotide (123G>C)
        boolean mutAA = token.matches(featureMutAA);
        boolean mutNT = token.matches(featureMutNT);
        if (mutAA || mutNT) {
            //We know it has the right form, but may
            //not be valid feature name or mutation
            //which exists.
            String[] items = token.toUpperCase().split(":");
            String name = items[0].trim().toUpperCase();
            String coords = items[1];
            int coordLength = coords.length();

            Map<Integer, BasicFeature> genomePosList;

            //Should never match both mutation notations
            if (mutAA) {
                String refSymbol = coords.substring(0, 1);
                String mutSymbol = coords.substring(coordLength - 1);

                String strLoc = coords.substring(1, coordLength - 1);
                int location = Integer.parseInt(strLoc) - 1;

                genomePosList = FeatureDB.getMutationAA(name, location + 1, refSymbol, mutSymbol, genome);
            } else if (mutNT) {
                //Exclude the "A>T" at end
                String strLoc = coords.substring(0, coordLength - 3);
                String refSymbol = coords.substring(coordLength - 3, coordLength - 2);
                int location = Integer.parseInt(strLoc) - 1;
                genomePosList = FeatureDB.getMutationNT(name, location + 1, refSymbol, genome);
            } else {
                //This should never happen
                throw new IllegalArgumentException("Something went wrong parsing input token");
            }

            for (int genomePos : genomePosList.keySet()) {
                Feature feature = genomePosList.get(genomePos);
                //Zoom in on mutation of interest
                //The +2 accounts for centering on the center of the amino acid, not beginning
                //and converting from 0-based to 1-based (which getStartEnd expects)
                int[] locs = getStartEnd("" + (genomePos + 2));
                return new SearchResult(ResultType.LOCUS, feature.getChr(), locs[0], locs[1]);
            }
        }

        return null;
    }

    private NamedFeature searchWebservice(String str) {
        try {
            List<Locus> results = SearchAPI.search(str, genome.getId());
            if (results != null && results.size() > 0) {

                results.sort((r1, r2) -> Integer.compare(r2.getLength(), r1.getLength()));

                // If the genome defines long chromosome names, prefer those
                Set<String> lnames = new HashSet<>(genome.getLongChromosomeNames());
                if (lnames.size() > 0) {
                    for (Locus l : results) {
                        if (lnames.contains(l.getChr())) {
                            return new BasicFeature(l.getChr(), l.getStart(), l.getEnd());
                        }
                    }
                }

                Locus l = results.get(0);
                String chr = genome == null ? l.getChr() : genome.getCanonicalChrName(l.getChr());
                return new BasicFeature(chr, l.getStart(), l.getEnd());
            }
        } catch (Exception e) {
            log.error("Search webservice error", e);
        }
        return null;
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
                Chromosome chromosome = genome.getChromosome(chr);
                if (chromosome == null) {
                    // try entire search string, chr name may have embedded colon
                    if (genome.getChromosome(searchString) != null) {
                        chr = searchString;
                        startEnd = null;
                    }
                } else {
                    String posString = searchString.substring(colonIdx).replace(":", "");
                    startEnd = getStartEnd(posString);

                }
            }
        }

        // Show the "All chromosomes" view if the search string is "*"
        if (chr.equals("*") || chr.toLowerCase().equals("all")) {
            return new SearchResult(ResultType.LOCUS, Globals.CHR_ALL, 0, Integer.MAX_VALUE);
        }

        //startEnd will have coordinates if found.
        Chromosome chromosome = genome.getChromosome(chr);
        //If we couldn't find chromosome, check
        //whole string
        if (chromosome == null) {
            chr = tokens[0];
            chromosome = genome.getChromosome(chr);
            if (chromosome != null) {
                //Found chromosome
                startEnd = null;
            }
        }

        if (chromosome != null && !searchString.equals(Globals.CHR_ALL)) {
            chr = chromosome.getName();
            if (startEnd == null) {
                return new SearchResult(ResultType.LOCUS, chr, 0, chromosome.getLength());
            } else {
                int start = Math.min(startEnd[0], startEnd[1]);
                int end = Math.max(startEnd[0], startEnd[1]);
                return new SearchResult(ResultType.LOCUS, chr, start, end);
            }
        }
        return null;
    }

    private void showFlankedRegion(String chr, int start, int end) {
        int flankingRegion = PreferencesManager.getPreferences().getAsInt(Constants.FLANKING_REGION);
        int delta;
        if ((end - start) == 1) {
            delta = 20; // Don't show flanking region for single base jumps, use 40bp window
        } else if (flankingRegion < 0) {
            delta = (-flankingRegion * (end - start)) / 100;
        } else {
            delta = flankingRegion;
        }
        start = Math.max(0, start - delta);
        end = end + delta;

        referenceFrame.jumpTo(chr, start, end);

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
            int start = Math.max(0, Integer.parseInt(startString) - 1);

            // Default value for end
            int end = start + 1;
            if (posTokens.length > 1) {
                String endString = posTokens[1].replaceAll(",", "");
                end = Integer.parseInt(endString);
            }

            if (posTokens.length == 1 || (end >= start && (end - start) < 10)) {
                int center = (start + end) / 2;
                int widen = 20;
                start = center - widen;
                start = Math.max(0, start);
                end = center + widen;
            }

            return new int[]{Math.min(start, end), Math.max(start, end)};
        } catch (NumberFormatException numberFormatException) {
            return null;
        }

    }

    public enum ResultType {
        FEATURE,
        LOCUS,
        ERROR,
        LIFTOVER
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

        public SearchResult() {
            this(ResultType.ERROR, null, -1, -1);
        }

        public SearchResult(ResultType type, String chr, int start, int end) {
            this.type = type;
            this.chr = chr;
            this.start = start;
            this.end = end;
            this.locus = Locus.getFormattedLocusString(chr, start, end);
        }

        public SearchResult(NamedFeature feature) {
            this(ResultType.FEATURE, feature.getChr(), feature.getStart(), feature.getEnd());
            this.feature = feature;
            this.locus = Locus.getFormattedLocusString(chr, start, end);
        }

        void setMessage(String message) {
            this.message = message;
        }

        public String getMessage() {
            return this.message;
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
                return this.locus;
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
                return feature.getName() + " (" + this.locus + ")";
            } else {
                return this.locus;
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

        //May be null
        @ForTesting
        public NamedFeature getFeature() {
            return feature;
        }
    }

    /**
     * Get a list of search results from the provided objects,
     * which must be IGVNamedFeature objects.
     *
     * @param objects
     * @return
     */
    public static List<SearchResult> getResults(List<IGVNamedFeature> objects) {
        List<SearchResult> results = new ArrayList<SearchResult>(objects.size());
        for (IGVNamedFeature f : objects) {
            results.add(new SearchCommand.SearchResult(f));
        }
        return results;
    }

    private static boolean isInteger(String str) {
        for (int i = 0; i < str.length(); i++) {
            char c = str.charAt(i);
            if (c < '0' || c > '9') return false;
        }
        return true;
    }
}