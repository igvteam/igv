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

package org.broad.igv.util;

import com.google.common.base.Function;
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


    @Test
    public void testGetIGVExtension() {

        String path = "/foo/bar/mydata.igv";
        assertEquals("igv", ParsingUtils.getIGVExtension(path));

        path = "/foo/bar/mydata.igv.gz";
        assertEquals("igv", ParsingUtils.getIGVExtension(path));

        path = "/foo/bar/mydata.igv.txt";
        assertEquals("igv", ParsingUtils.getIGVExtension(path));

        path = "/foo/bar/mydata.igv.xls";
        assertEquals("igv", ParsingUtils.getIGVExtension(path));

        path = "/foo/bar/mydata.igv.txt.gz";
        assertEquals("igv", ParsingUtils.getIGVExtension(path));

    }

    //@Test
    public void compareSpeedPatternDotSplit() throws Exception {
        int nTrials = 500000;
        TestStringSupplier supplier = new TestStringSupplier(nTrials);

        final char cdelim = '\t';
        final String sdelim = String.valueOf(cdelim);

        Function<String, Void> patternSplitPredicate = new Function<String, Void>() {
            @Override
            public Void apply(String input) {
                String[] tokens = Globals.tabPattern.split(input);
                return null;
            }
        };

        //ParsingUtils.split seems to be about 2x as fast, that is
        //takes 1/2 the time
        Function<String, Void> parsingUtilsPredicate = new Function<String, Void>() {
            String[] buffer = new String[20];

            @Override
            public Void apply(String input) {
                int count = htsjdk.tribble.util.ParsingUtils.split(input, buffer, cdelim);
                return null;
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

        Function<String[], Void> stringBuilderFunc = new Function<String[], Void>() {
            @Override
            public Void apply(String[] input) {
                StringBuilder stringBuilder = new StringBuilder();
                for (int el = 0; el < input.length; el++) {
                    stringBuilder.append(input[el]);
                    if (el < input.length - 1) {
                        stringBuilder.append(sdelim);
                    }
                }
                return null;
            }
        };

        //ParsingUtils.split seems to be about 2x as fast, that is
        //takes 1/2 the time
        Function<String[], Void> stringAddFunc = new Function<String[], Void>() {
            @Override
            public Void apply(String[] input) {
                String result = "";

                for (int el = 0; el < input.length; el++) {
                    result += input[el];
                    if (el < input.length - 1) {
                        result += sdelim;
                    }
                }

                return null;

                //String res = htsjdk.tribble.util.ParsingUtils.join(sdelim, input);
                //return true;
            }
        };

        supplier.reset();
        System.out.println("\nStringBuilder");
        TestUtils.timeMethod(supplier, stringBuilderFunc, nTrials);

        supplier.reset();
        System.out.println("\nStringAdd");
        TestUtils.timeMethod(supplier, stringAddFunc, nTrials);
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

