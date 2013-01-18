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

package org.broad.igv.util;

import com.google.common.base.Predicate;
import com.google.common.base.Supplier;
import junit.framework.Assert;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.Globals;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.WindowFunction;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import java.awt.*;
import java.util.Calendar;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

/**
 * User: jrobinso
 * Date: Feb 8, 2010
 */
public class ParsingUtilsTest extends AbstractHeadlessTest {

    public final static String characters = "0123456789abcdefghijklmnopqrstuvwxyz";
    public final static int numChars = characters.length();

    private static final int connectTimeout = 30 * 1000;
    @Rule
    public TestRule testTimeout = new Timeout(4 * connectTimeout);

    @Before
    public void setUp() throws Exception {
        super.setUp();
        Globals.CONNECT_TIMEOUT = connectTimeout;
    }

    private String genRandString() {
        int numWords = 10;
        int max_length = 20;
        String ret = "";
        for (int _ = 0; _ < numWords; _++) {
            ret += getRandWord(max_length) + "\t";
        }
        return ret;
    }

    private String getRandWord(int max_length) {
        int length = (int) Math.random() * max_length + 1;
        String ret = "";
        for (int _ = 0; _ < length; _++) {
            ret += characters.charAt((int) Math.random() * numChars);
        }
        return ret;
    }

    @Test
    public void testSplit1() {
        String blankColumnLine = "a\tb\t\td";
        String[] tokens = Globals.tabPattern.split(blankColumnLine);
        int nTokens = tokens.length;
        assertEquals(4, nTokens);
        assertEquals("a", tokens[0]);
        assertEquals("b", tokens[1]);
        assertEquals("", tokens[2]);
        assertEquals("d", tokens[3]);
    }

    @Test
    public void testSplit2() {
        String blankColumnLine = "a\tb\t\td\t";
        String[] tokens = Globals.tabPattern.split(blankColumnLine);
        int nTokens = tokens.length;
        assertEquals(4, nTokens);
        assertEquals("a", tokens[0]);
        assertEquals("b", tokens[1]);
        assertEquals("", tokens[2]);
        assertEquals("d", tokens[3]);
    }

    @Test
    public void testSplit3() {
        String blankColumnLine = "\ta\t\tb\t\t\td\t";
        String[] tokens = Globals.tabPattern.split(blankColumnLine);
        int nTokens = tokens.length;
        String[] expTokens = new String[]{"", "a", "", "b", "", "", "d"};
        assertEquals(expTokens.length, nTokens);
        int ii = 0;
        for (String exp : expTokens) {
            Assert.assertEquals(exp, tokens[ii]);
            ii++;
        }
    }

    @Test
    public void testGetContentLengthFTP() {
        long contLength = ParsingUtils.getContentLength(TestUtils.AVAILABLE_FTP_URL);
        assertTrue("Error retrieving content length: " + contLength, contLength > 0);

        long start_time = System.currentTimeMillis();
        assertTrue(ParsingUtils.getContentLength(TestUtils.UNAVAILABLE_FTP_URL) == -1);
        long end_time = System.currentTimeMillis();
        assertTrue(end_time - start_time < Globals.CONNECT_TIMEOUT + 1000);
    }


    @Test
    public void testGetLastModified_HTTP() throws Exception{
        tstGetLastModified(HttpUtilsTest.broadURLString);
    }

    @Test
    public void testGetLastModified_FTP() throws Exception{
        tstGetLastModified(TestUtils.AVAILABLE_FTP_URL);
    }

    @Test
    public void testGetLastModified_File() throws Exception{
        tstGetLastModified(TestUtils.DATA_DIR + "bed/test.bed");
    }

    private void tstGetLastModified(String path) {
        long modD = ParsingUtils.getLastModified(path);
        assertTrue(modD > 0);

        //Assuming causality is still intact
        assertTrue(modD < Calendar.getInstance().getTime().getTime());
    }

    @Test
    public void testParseInt() {
        String with_commas = "123456";
        int expected = 123456;
        int actual = ParsingUtils.parseInt(with_commas);
        assertEquals(expected, actual);

        String exp_not = "3.5e4";
        expected = 35000;
        assertEquals(expected, ParsingUtils.parseInt(exp_not));
    }

    @Test
    public void testParseTrackLine() {
        String trackLine = "track type=bigWig name=\"Track 196\" visibility=2 " +
                "description=\" CD34 - H3K27me3 - hg19 - 18.7 M/20.9 M - 61P7DAAXX.6\" " +
                "maxHeightPixels=70 viewLimits=0:18 windowingFunction=mean autoScale=off " +
                "bigDataUrl=http://www.broadinstitute.org/epigenomics/dataportal/track_00196.portal.bw " +
                "color=255,0,0";

        TrackProperties props = new TrackProperties();
        ParsingUtils.parseTrackLine(trackLine, props);
        assertEquals("Track 196", props.getName());
        assertEquals(Track.DisplayMode.EXPANDED, props.getDisplayMode());
        assertEquals(" CD34 - H3K27me3 - hg19 - 18.7 M/20.9 M - 61P7DAAXX.6", props.getDescription());
        assertEquals(70, props.getHeight());
        assertEquals(0, props.getMinValue(), 1.0e-9);
        assertEquals(18, props.getMaxValue(), 1.0e-9);
        assertEquals(WindowFunction.mean, props.getWindowingFunction());
        assertEquals(false, props.isAutoScale());
        assertEquals(new Color(255, 0, 0), props.getColor());
        assertEquals("http://www.broadinstitute.org/epigenomics/dataportal/track_00196.portal.bw", props.getDataURL());
    }

    //@Test
    public void compareSpeedPatternDotSplit() throws Exception {
        int nTrials = 500000;
        TestStringSupplier supplier = new TestStringSupplier(nTrials);

        final char cdelim = '\t';
        final String sdelim = String.valueOf(cdelim);

        Predicate<String> patternSplitPredicate = new Predicate<String>() {
            @Override
            public boolean apply(String input) {
                String[] tokens = Globals.tabPattern.split(input);
                return true;
            }
        };

        //ParsingUtils.split seems to be about 2x as fast, that is
        //takes 1/2 the time
        Predicate<String> parsingUtilsPredicate = new Predicate<String>() {
            String[] buffer = new String[20];

            @Override
            public boolean apply(String input) {
                int count = org.broad.tribble.util.ParsingUtils.split(input, buffer, cdelim);
                return true;
            }
        };

        supplier.reset();
        System.out.println("\nPattern.split");
        TestUtils.timeMethod(supplier, patternSplitPredicate, nTrials);

        supplier.reset();
        System.out.println("\nParsingUtils.split");
        TestUtils.timeMethod(supplier, parsingUtilsPredicate, nTrials);
    }

    //@Test
    public void compareSpeedStringJoin() throws Exception {
        int nTrials = 5000;
        TestStringArraySupplier supplier = new TestStringArraySupplier(nTrials, 100);

        final char cdelim = '\t';
        final String sdelim = String.valueOf(cdelim);

        Predicate<String[]> stringBuilderPredicate = new Predicate<String[]>() {
            @Override
            public boolean apply(String[] input) {
                StringBuilder stringBuilder = new StringBuilder();
                for (int el = 0; el < input.length; el++) {
                    stringBuilder.append(input[el]);
                    if (el < input.length - 1) {
                        stringBuilder.append(sdelim);
                    }
                }
                return true;
            }
        };

        //ParsingUtils.split seems to be about 2x as fast, that is
        //takes 1/2 the time
        Predicate<String[]> stringAddPred = new Predicate<String[]>() {
            @Override
            public boolean apply(String[] input) {
                String result = "";

                for (int el = 0; el < input.length; el++) {
                    result += input[el];
                    if (el < input.length - 1) {
                        result += sdelim;
                    }
                }

                return true;

                //String res = org.broad.tribble.util.ParsingUtils.join(sdelim, input);
                //return true;
            }
        };

        supplier.reset();
        System.out.println("\nStringBuilder");
        TestUtils.timeMethod(supplier, stringBuilderPredicate, nTrials);

        supplier.reset();
        System.out.println("\nStringAdd");
        TestUtils.timeMethod(supplier, stringAddPred, nTrials);
    }

    private class TestStringArraySupplier implements Supplier<String[]> {
        final String[][] testStringArrays;
        private int counter = 0;

        TestStringArraySupplier(int nTrials, int maxStringsPerArray) {
            testStringArrays = new String[nTrials][];
            for (int ii = 0; ii < nTrials; ii++) {
                int numEntries = (int) Math.floor(Math.random() * maxStringsPerArray);
                testStringArrays[ii] = genTestStringArray(numEntries);
            }
        }

        @Override
        public String[] get() {
            return testStringArrays[counter++];
        }

        public void reset() {
            counter = 0;
        }
    }

    private class TestStringSupplier implements Supplier<String> {

        final String[] testStrings;
        private int counter = 0;

        TestStringSupplier(int nTrials) {
            testStrings = genTestStringArray(nTrials);
        }


        @Override
        public String get() {
            return testStrings[counter++];
        }

        public void reset() {
            counter = 0;
        }

    }

    private String[] genTestStringArray(int numEntries) {
        String[] testStrings = new String[numEntries];

        //Generate test data
        for (int ii = 0; ii < numEntries; ii++) {
            testStrings[ii] = genRandString();
        }
        return testStrings;
    }


}

