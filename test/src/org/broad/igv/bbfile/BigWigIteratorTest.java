/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
        int nTrials = 10;
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
