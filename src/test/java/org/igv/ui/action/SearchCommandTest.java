package org.igv.ui.action;

import htsjdk.tribble.NamedFeature;
import junit.framework.AssertionFailedError;
import org.igv.AbstractHeadlessTest;
import org.igv.feature.BasicFeature;
import org.igv.feature.genome.ChromAlias;
import org.igv.feature.genome.ChromAliasSource;
import org.igv.feature.genome.Genome;
import org.igv.util.TestUtils;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;
import java.util.*;

import static junit.framework.Assert.assertEquals;
import static org.junit.Assert.*;


/**
 * User: jacob
 * Date: 2011/12/19
 */
public class SearchCommandTest extends AbstractHeadlessTest {

    @Before
    public void setUp() throws Exception {
        super.setUp();
    }

    @Test
    public void testSingleLOCUSs() throws Exception {

        int min = 1;
        int max = 22;
        String[] chrs = new String[max - min + 1];
        String[] nums = new String[max - min + 1];
        String chr;

        for (int cn = min; cn <= max; cn++) {
            chr = "chr" + cn;
            chrs[cn - min] = chr;
            nums[cn - min] = "" + cn;
        }
        tstFeatureTypes(chrs, SearchCommand.ResultType.LOCUS);
        tstFeatureTypes(nums, SearchCommand.ResultType.LOCUS);
    }

    @Test
    public void testChromoWithColon() throws Exception {


        String [] chrNames = {"chr10", "chr1", "chr20"};
        String [][] aliases = {{"abc:123", "abc:456", "abc:789"}, {"xy:12"}, {"aaa:bbb"}};
        List<List<String>> synonymsList = new ArrayList<>();
        synonymsList.add(Arrays.asList("chr10", "abc:123", "abc:456", "abc:789"));
        synonymsList.add(Arrays.asList("chr1", "xy:12"));
        synonymsList.add(Arrays.asList("chr20", "aaa:bbb"));
        genome.addChrAliases(synonymsList);


        SearchCommand cmd;
        for (int i = 0; i < chrNames.length; i++) {

            String chr = chrNames[i];
            String [] synonyms = aliases[i];

            for (String searchStr : synonyms) {
                cmd = new SearchCommand(null, searchStr, genome);
                List<SearchCommand.SearchResult> results = cmd.runSearch(cmd.searchString);
                assertEquals(1, results.size());
                assertEquals(SearchCommand.ResultType.LOCUS, results.get(0).type);
                assertEquals(chr, results.get(0).chr);
            }
        }

        String[] queries = new String[]{"X:100-1000", "Y:100-1000", "Y:12", "1\t 100", "X\t50", "X\t50\t500"};
        tstFeatureTypes(queries, SearchCommand.ResultType.LOCUS);

        genome = TestUtils.loadGenome();
    }

    @Test
    public void testSingleFeatures() throws Exception {
        String[] features = {"EGFR", "ABO", "BRCA1", "egFr"};
        tstFeatureTypes(features, SearchCommand.ResultType.FEATURE);
    }

    /**
     * Test that mutations are parsed correctly. They get returned
     * as type LOCUS.
     *
     * @throws Exception
     */
    @Test
    @Ignore("Fails unless tests are run in separate JVMs")
    public void testFeatureMuts() throws Exception {
        String[] features = {"EGFR:M1I", "EGFR:G5R", "egfr:g5r", "egfr:r2*"};
        tstFeatureTypes(features, SearchCommand.ResultType.LOCUS);

        String egfr_seq = "atgcgaccctccgggacggccggggcagcgctcctggcgctgctggctgcgctc";
        String[] targets = new String[]{"A", "G", "C", "T", "a", "g", "c", "t"};
        String mut_let, tar_let;
        for (int let = 0; let < egfr_seq.length(); let++) {
            mut_let = egfr_seq.substring(let, let + 1);
            if (Math.random() > 0.5) {
                mut_let = mut_let.toUpperCase();
            }
            tar_let = targets[let % targets.length];
            String[] features2 = {"EGFR:" + (let + 1) + mut_let + ">" + tar_let};
            tstFeatureTypes(features2, SearchCommand.ResultType.LOCUS);
        }
    }

    /*
    Run a bunch of queries, test that they all return the same type
     */
    private void tstFeatureTypes(String[] queries, SearchCommand.ResultType type) {
        SearchCommand cmd;
        for (String f : queries) {
            cmd = new SearchCommand(null, f, genome);
            List<SearchCommand.SearchResult> results = cmd.runSearch(cmd.searchString);
            try {
                assertTrue(containsType(results, type));
            } catch (AssertionFailedError e) {
                System.out.println(f + " :" + results.get(0).getMessage());
                throw e;
            }
        }
    }

    private boolean containsType(List<SearchCommand.SearchResult> results, SearchCommand.ResultType check) {
        boolean contains = false;
        for (SearchCommand.SearchResult result : results) {
            contains |= result.type == check;
        }
        return contains;
    }

    @Test
    public void testMultiFeatures() throws Exception {
        tstMultiFeatures(" ");
        tstMultiFeatures("  ");
        tstMultiFeatures("   ");
        tstMultiFeatures("    ");
        tstMultiFeatures("\t");
        tstMultiFeatures("\r\n"); //This should never happen
    }

    public void tstMultiFeatures(String delim) throws Exception {

        //Not sure if this is correct behavior or not
        String[] tokens = {"EgfR", "ABO", "BRCA1", "chr1:1-100"};
        SearchCommand.ResultType[] types = new SearchCommand.ResultType[]{
                SearchCommand.ResultType.FEATURE,
                SearchCommand.ResultType.FEATURE,
                SearchCommand.ResultType.FEATURE,
                SearchCommand.ResultType.LOCUS,
                SearchCommand.ResultType.LOCUS,
                SearchCommand.ResultType.LOCUS
        };
        String searchStr = tokens[0];
        for (int ii = 1; ii < tokens.length; ii++) {
            searchStr += delim + tokens[ii];
        }
        SearchCommand cmd;
        cmd = new SearchCommand(null, searchStr, genome);
        List<SearchCommand.SearchResult> results = cmd.runSearch(cmd.searchString);

        for (int ii = 0; ii < tokens.length; ii++) {
            SearchCommand.SearchResult result = results.get(ii);
            try {
                assertEquals(types[ii], result.type);
                assertEquals(tokens[ii].toLowerCase(), result.getShortName().toLowerCase());
            } catch (AssertionFailedError e) {
                System.out.println(searchStr + " :" + result.getMessage());
                throw e;
            }
        }
    }

    @Test
    public void testMultiLOCUSs() throws Exception {
        String[] tokens = {"chr1", "chr5", "4", "12", "X", "Y"};
        String searchStr = tokens[0];
        for (int ii = 1; ii < tokens.length; ii++) {
            searchStr += " " + tokens[ii];
        }

        SearchCommand cmd;
        cmd = new SearchCommand(null, searchStr, genome);
        List<SearchCommand.SearchResult> results = cmd.runSearch(cmd.searchString);

        for (int ii = 0; ii < tokens.length; ii++) {
            SearchCommand.SearchResult result = results.get(ii);

            assertEquals(SearchCommand.ResultType.LOCUS, result.type);
            assertTrue(result.getLocus().contains(tokens[ii]));
        }
    }



    /**
     * Saving because we might want these test cases, but don't
     * want to test checkTokenType specifically
     */
    @Deprecated
    public void testTokenChecking() {
        String[] chromos = {"chr3", "chr20", "chrX", "chrY"};
        SearchCommand cmd = new SearchCommand(null, "", genome);
        for (String chr : chromos) {
            assertEquals(SearchCommand.ResultType.LOCUS, cmd.checkTokenType(chr));
        }

        String[] starts = {"39,239,480", "958392", "0,4829,44", "5"};
        String[] ends = {"40,321,111", "5", "48153181,813156", ""};
        for (int ii = 0; ii < starts.length; ii++) {
            String tstr = chromos[ii] + ":" + starts[ii] + "-" + ends[ii];
            assertEquals(SearchCommand.ResultType.LOCUS, cmd.checkTokenType(tstr));
            tstr = chromos[ii] + "\t" + starts[ii] + "  " + ends[ii];
            assertEquals(SearchCommand.ResultType.LOCUS, cmd.checkTokenType(tstr));
        }

        String[] errors = {"egfr:1-100", "   ", "chr1\t1\t100\tchr2"};
        for (String s : errors) {
            //System.out.println(s);
            assertEquals(SearchCommand.ResultType.ERROR, cmd.checkTokenType(s));
        }
    }


    @Test
    public void testMutationSearch() throws Exception {

        String name = "EGFR";
        // EGFR starts with proteins MRPSG
        String[] symbols = new String[]{"M", "R", "P", "S", "G"};
        //All of these should be possible with a SNP from the EGFR sequence
        String[] muts = new String[]{"I", "R", "P", "P", "R"};
        Map<Integer, BasicFeature> matches;
        List<NamedFeature> possibles = genome.getFeatureDB().getFeaturesMatching(name);
        for (int ii = 0; ii < symbols.length; ii++) {
            matches = SearchCommand.getMutationAA(possibles, ii + 1, symbols[ii], muts[ii], genome);
            assertEquals(1, matches.size());
            for (int pos : matches.keySet()) {
                assertEquals(name, matches.get(pos).getName());
            }
        }

        name = "EGFLAM";
        int exp_start = 38439399;
        possibles = genome.getFeatureDB().getFeaturesMatching(name);
        matches = SearchCommand.getMutationAA(possibles, 2, "H", "H", genome);
        assertEquals(1, matches.size());
        for (int geneloc : matches.keySet()) {
            assertEquals(exp_start, matches.get(geneloc).getStart());
        }

        String[] others = new String[]{"I", "M", "T"};
        for (String c : others) {
            matches = SearchCommand.getMutationAA(possibles, 2, "H", c, genome);
            assertEquals(0, matches.size());
        }
    }

    @Test
    public void testMutationSearchNegStrand() throws Exception {
        String name = "KRAS";
        int exp_start = 25249446;
        List<NamedFeature> possibles = genome.getFeatureDB().getFeaturesMatching(name);
        Map<Integer, BasicFeature> matches =SearchCommand.getMutationAA(possibles, 1, "M", "I", genome);
        assertEquals(1, matches.size());
        for (int geneloc : matches.keySet()) {
            assertEquals(exp_start, matches.get(geneloc).getStart());
        }

    }

    @Test
    public void testMutationSearchFail() throws Exception {
        String name = "EGFR";
        String[] symbols = "R,P,S,G,M".split(",");
        List<NamedFeature> possibles = genome.getFeatureDB().getFeaturesMatching(name);
        for (int ii = 0; ii < symbols.length; ii++) {
            Map<Integer, BasicFeature>  matches = SearchCommand.getMutationAA(possibles, ii + 1, symbols[ii], "M", genome);
            assertEquals(0, matches.size());
        }
    }

    @Test
    public void testMutationSearchNT() throws Exception {
        String name = "EGFR";
        String[] bps = new String[]{"A", "T", "G"};
        List<NamedFeature> possibles = genome.getFeatureDB().getFeaturesMatching(name);
        for (int ii = 0; ii < bps.length; ii++) {
            Map<Integer, BasicFeature> matches = SearchCommand.getMutationNT(possibles, ii + 1, bps[ii], genome);
            assertEquals(1, matches.size());
        }
    }

    @Test
    public void testMutationSearchNTNegStrand() throws Exception {
        String name = "KRAS";
        String[] bps = new String[]{"A", "T", "G"};
        List<NamedFeature> possibles = genome.getFeatureDB().getFeaturesMatching(name);

        for (int ii = 0; ii < bps.length; ii++) {
            Map<Integer, BasicFeature> matches = SearchCommand.getMutationNT(possibles, ii + 1, bps[ii], genome);
            assertEquals(1, matches.size());
        }

        //Exon 3 of KRAS, starting at amino acid 56
        int startNT = (56 - 1) * 3;
        char[] bps2 = "CTCGACACAGCAGGT".toCharArray();
        for (int ii = 0; ii < bps2.length; ii++) {
            Map<Integer, BasicFeature> matches = SearchCommand.getMutationNT(possibles, ii + startNT + 1, String.valueOf(bps2[ii]), genome);
            assertEquals(1, matches.size());
        }
    }
}


