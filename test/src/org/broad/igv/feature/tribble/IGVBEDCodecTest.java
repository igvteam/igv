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

package org.broad.igv.feature.tribble;

import com.google.common.base.Predicate;
import com.google.common.base.Supplier;
import org.apache.commons.lang.StringUtils;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.Globals;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012-Dec-11
 */
public class IGVBEDCodecTest extends AbstractHeadlessTest {


    @Test
    public void decodeSpeedTest() throws Exception {
        int nTrials = 5000000;
        BEDStringSupplier supplier = new BEDStringSupplier(nTrials);
        final IGVBEDCodec codec = new IGVBEDCodec();

        Predicate<String> decodePredicate = new Predicate<String>() {
            @Override
            public boolean apply(String input) {
                BasicFeature feat = codec.decode(input);
                return true;
            }
        };

        supplier.reset();
        String jVersion = System.getProperty(Globals.JAVA_VERSION_STRING);
        System.out.println("\nIGVBEDCodec.decode. java version " + jVersion);
        long[] times = TestUtils.timeMethod(supplier, decodePredicate, nTrials);
        long median = times[times.length / 2];

        long benchTime = TestUtils.getBenchmarkTime();

        int maxMultiplier = 200000;
        if (jVersion.contains("1.7")) {
            maxMultiplier = 10000;
        }
        assertTrue("Decoding median speed too slow", median < benchTime / maxMultiplier);
    }

    private class BEDStringSupplier implements Supplier<String> {

        //final String[] testStrings;
        private int counter = 0;


        private final String chr = "chr1";
        private final int maxLength = 500;

        BEDStringSupplier(int nTrials) {
//            testStrings = new String[nTrials];
//            for(int ft=0; ft < nTrials; ft++){
//                testStrings[ft] = gen(ft);
//            }
        }

        private String gen(int start) {
            String end = "" + (start + (int) Math.random() * maxLength);
            String[] dat = {chr, "" + start, end, "0", "0", "+"};
            return StringUtils.join(dat, '\t');
        }


        @Override
        public String get() {
            return gen(counter++);
            //return testStrings[counter++];
        }

        public void reset() {
            counter = 0;
        }

    }
}
