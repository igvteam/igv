/*
 * Copyright (c) 2007-2012 by The Broad Institute of MIT and Harvard.All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.ui.action;

import junit.framework.AssertionFailedError;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.genome.GenomeImpl;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.*;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


/**
 * User: jacob
 * Date: 2011/12/19
 */
public class SearchCommandTest extends AbstractHeadlessTest {

    @Test
    public void testSingleChromosomes() throws Exception {

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
        tstFeatureTypes(chrs, SearchCommand.ResultType.CHROMOSOME);
        tstFeatureTypes(nums, SearchCommand.ResultType.CHROMOSOME);
    }

    @Test
    public void testChromoWithColon() throws Exception {

        List<String> chrNames = Arrays.asList("chr10", "chr1", "chr20");
        List<List<String>> aliases = Arrays.asList(
                Arrays.asList("abc:123", "abc:456", "abc:789"),
                Arrays.asList("xy:12"),
                Arrays.asList("aaa:bbb"));


        Collection<Collection<String>> synonymsList = new ArrayList<Collection<String>>();
        for (int i = 0; i < chrNames.size(); i++) {
            List<String> synonyms = new ArrayList<String>();
            synonyms.addAll(aliases.get(i));
            synonyms.add(chrNames.get(i));
            synonymsList.add(synonyms);
        }

        if (genome instanceof GenomeImpl) {
            ((GenomeImpl) genome).addChrAliases(synonymsList);
        }

        SearchCommand cmd;
        for (int i = 0; i < chrNames.size(); i++) {

            String chr = chrNames.get(i);
            List<String> synonyms = aliases.get(i);

            for (String searchStr : synonyms) {
                cmd = new SearchCommand(null, searchStr, genome);
                List<SearchCommand.SearchResult> results = cmd.runSearch(cmd.searchString);
                assertEquals(1, results.size());
                assertEquals(SearchCommand.ResultType.CHROMOSOME, results.get(0).type);
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
                SearchCommand.ResultType.CHROMOSOME,
                SearchCommand.ResultType.CHROMOSOME
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
                assertEquals(tokens[ii].toLowerCase(), result.getLocus().toLowerCase());
            } catch (AssertionFailedError e) {
                System.out.println(searchStr + " :" + result.getMessage());
                throw e;
            }
        }
    }

    @Test
    public void testMultiChromosomes() throws Exception {
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

            assertEquals(SearchCommand.ResultType.CHROMOSOME, result.type);
            assertTrue(result.getLocus().contains(tokens[ii]));
        }
    }

    @Test
    public void testError() throws Exception {
        String[] tokens = {"ueth", "EGFRa", "BRCA56", "EGFR:?1?"};
        tstFeatureTypes(tokens, SearchCommand.ResultType.ERROR);
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
            assertEquals(SearchCommand.ResultType.CHROMOSOME, cmd.checkTokenType(chr));
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

}
