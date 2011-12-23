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
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tools.IgvTools;
import org.junit.Before;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


/**
 * User: jacob
 * Date: 2011/12/19
 */
public class SearchCommandTest {

    private Genome genome;

    @Before
    public void setUp() throws Exception {
        String dataFileName = "test/data/genomes/hg18.genome";
        Globals.setHeadless(true);
        this.genome = IgvTools.loadGenome(dataFileName, true);
    }

    @Test
    public void testSingleChromosomes() throws Exception {

        int min = 1;
        int max = 22;
        String[] chrs = new String[max - min + 1];
        String chr;

        for (int cn = min; cn <= max; cn++) {
            chr = "chr" + cn;
            chrs[cn - min] = chr;
        }
        tstFeatureTypes(chrs, SearchCommand.ResultType.CHROMOSOME);
    }

    @Test
    public void testSingleFeatures() throws Exception {
        String[] features = {"EGFR", "ABO", "BRCA1"};
        tstFeatureTypes(features, SearchCommand.ResultType.FEATURE);
    }

    /*
    Run a bunch of queries, test that they all return the same type
     */
    private void tstFeatureTypes(String[] queries, SearchCommand.ResultType type) {
        SearchCommand cmd;
        for (String f : queries) {
            cmd = new SearchCommand(null, f, genome);
            SearchCommand.SearchResult result = cmd.runSearch(cmd.searchString).get(0);
            try {
                assertEquals(type, result.type);
            } catch (AssertionFailedError e) {
                System.out.println(f + " :" + result.getMessage());
                throw e;
            }
        }
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
        String[] tokens = {"EgfR", "ABO", "BRCA1", "chr1:1-100"};
        SearchCommand.ResultType[] types = new SearchCommand.ResultType[]{
                SearchCommand.ResultType.FEATURE,
                SearchCommand.ResultType.FEATURE,
                SearchCommand.ResultType.FEATURE,
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
                assertEquals(tokens[ii].toLowerCase(), result.getLocus().toLowerCase());
            } catch (AssertionFailedError e) {
                System.out.println(searchStr + " :" + result.getMessage());
                throw e;
            }
        }
    }

    @Test
    public void testMultiChromosomes() throws Exception {
        String[] tokens = {"chr1", "chr5"};
        String searchStr = tokens[0];
        for (int ii = 1; ii < tokens.length; ii++) {
            searchStr += " " + tokens[ii];
        }

        SearchCommand cmd;
        cmd = new SearchCommand(null, searchStr, genome);
        List<SearchCommand.SearchResult> results = cmd.runSearch(cmd.searchString);

        for (int ii = 0; ii < tokens.length; ii++) {
            SearchCommand.SearchResult result = results.get(ii);
            //We add coordinates, and expect that loci may be
            //chr1:1-xxx where xxx is some large number.
            assertEquals(SearchCommand.ResultType.LOCUS, result.type);
            assertTrue(result.getLocus().contains(tokens[ii]));
        }
    }

    @Test
    public void testError() throws Exception {
        String[] tokens = {"ueth", "EGFRa", "BRCA56"};
        tstFeatureTypes(tokens, SearchCommand.ResultType.ERROR);
    }

/*    @Test
    public void testTokenParsing(){
        String[] chromos = {"chr3", "chr20","chrX", "chrY"};
        SearchCommand cmd = new SearchCommand(null, "", genome);
        for(String chr: chromos){
            assertEquals(SearchCommand.TokenType.CHROMOSOME, cmd.parseIndividualToken(chr));
        }
        
        String[] starts = {"39,239,480", "958392", "0,4829,44"};
        String[] ends = {"40,321,111","5","48153181,813156"};
        for(int ii=0; ii < starts.length; ii++){
            String tstr = chromos[ii] + ":" + starts[ii]+ "-" + ends[ii];
            assertEquals(SearchCommand.TokenType.LOCUS, cmd.parseIndividualToken(tstr));
            tstr = chromos[ii] + "\t" + starts[ii]+ "  " + ends[ii];
            assertEquals(SearchCommand.TokenType.LOCUS, cmd.parseIndividualToken(tstr));
        }
        
        String[] errors = {"egfr:1-100", "chr100", "chrz", "chr1:3","   "};
        for(String s: errors){
            System.out.println(s);
            assertEquals(SearchCommand.TokenType.ERROR, cmd.parseIndividualToken(s));
        }
    }*/

}
