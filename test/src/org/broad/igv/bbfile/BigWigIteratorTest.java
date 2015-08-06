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

package org.broad.igv.bbfile;

import com.google.common.base.Function;
import com.google.common.base.Supplier;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.Globals;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

/**
 * @author jacob
 * @date 2013-Jul-09
 */
public class BigWigIteratorTest extends AbstractHeadlessTest{

    @Test
    public void iterateSpeedTestLocal() throws Exception{
        String path = TestUtils.DATA_DIR + "wig/dummy_var_sample.bigwig";
        int nTrials = 500;
        iterateSpeedTst(path, nTrials, 2, 5);
    }

    @Test
    public void iterateSpeedTestRemote() throws Exception{
        String path = "ftp://ftp.broadinstitute.org/pub/igv/TEST/dummy_var_sample.bigwig";
        int nTrials = 50;
        iterateSpeedTst(path, nTrials, 2, 5);
    }

    public void iterateSpeedTst(String path, int nTrials, int j6Mult, int j7Mult) throws Exception{

        final BBFileReader reader = new BBFileReader(path);

        assertTrue(reader.isBigWigFile());

        long benchTime = TestUtils.getBenchmarkTime();


        Supplier<String> supplier = new Supplier<String>() {
            @Override
            public String get() {
                return "chr8";
            }
        };

        final Function<String, Void> func = new Function<String, Void>() {
            @Override
            public Void apply(String chr) {
                BigWigIterator iter = reader.getBigWigIterator(chr, 0, chr, Integer.MAX_VALUE, true);
                while(iter.hasNext()){
                    iter.next();
                }
                return null;
            }
        };

//        Predicate<String> pred = new Predicate<String>() {
//            @Override
//            public boolean apply(String input) {
//                func.apply(input);
//                return true;
//            }
//        };

        String jVersion = System.getProperty(Globals.JAVA_VERSION_STRING);
        System.out.println("\nBigWigIterator. java version " + jVersion);

        long[] times = TestUtils.timeMethod(supplier, func, nTrials);
        //Calculate average (in nanoseconds)
        double average = TestUtils.average(times);
        long median = times[times.length / 2];

        int maxMultiplier = j6Mult;
        if (jVersion.contains("1.7")) {
            maxMultiplier = j7Mult;
        }
        //we are somewhat forgiving, for the sake of portability.
        assertTrue("Median speed too slow", median < (benchTime / maxMultiplier) || median < 1e7);
    }
}
