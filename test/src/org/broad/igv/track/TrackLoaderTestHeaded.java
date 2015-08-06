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

package org.broad.igv.track;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.AbstractHeadedTest;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 * @author Jim Robinson
 * @date 10/3/11
 */
public class TrackLoaderTestHeaded extends AbstractHeadedTest {

    TrackLoader trackLoader;

    @Before
    public void setUp() throws Exception {
        super.setUp();
        trackLoader = new TrackLoader();
    }

    @Test
    public void testFilesHeaded() throws Exception {

        String ex_filename = "/vcf/example4-last-gsnap-2_fixed.vcf";
        Genome genome = TestUtils.loadGenome();
        List<String> finames = new ArrayList<String>(Arrays.asList(TrackLoaderTest.filenamesTryIndex));

        finames.add(ex_filename);

        for (String finame : finames) {
            TrackLoaderTest.tstLoadFi(trackLoader, TestUtils.DATA_DIR + finame, null, genome, true);
        }
    }


    @Test
    public void testReadVCF() throws Exception {
        String file = TestUtils.DATA_DIR + "vcf/outputPileup.flt1.vcf";
        TestUtils.createIndex(file);
        ResourceLocator locator = new ResourceLocator(file);
        //For files with 1 record, this threw a null pointer exception prior to r1595
        List<Track> tracks = igv.load(locator);

        Assert.assertEquals(1, tracks.size());
    }

    @Test
    public void testLoadBAMFtpHeaded() throws Exception{
        String path = "ftp://ftp.broadinstitute.org/pub/igv/TEST/HG00171.hg18.bam";
        Integer expected_tracks = 2;

        TrackLoaderTest.tstLoadFi(trackLoader, path, expected_tracks, genome, false);
    }
}
