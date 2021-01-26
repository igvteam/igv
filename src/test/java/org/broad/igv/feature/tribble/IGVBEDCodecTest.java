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

package org.broad.igv.feature.tribble;

import com.google.common.base.Function;
import com.google.common.base.Supplier;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import junit.framework.Assert;
import org.apache.commons.lang.StringUtils;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.Globals;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Ignore;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012-Dec-11
 */

public class IGVBEDCodecTest extends AbstractHeadlessTest {

    @Test
    public void testTabDelimited() throws Exception {
        String line = "chr1\t1051161\t1051177\tstrong GGGCGGGTGGGGCGGG\t0.81";
        IGVBEDCodec codec = new IGVBEDCodec();
        IGVFeature feature =  codec.decode(line);
        assertEquals("strong GGGCGGGTGGGGCGGG", feature.getName());
    }

    @Test
    public void testMultiTabDelimited() throws Exception {
        String line = "chr1\t1051161\t\t\t1051177\tstrong GGGCGGGTGGGGCGGG\t0.81";
        IGVBEDCodec codec = new IGVBEDCodec();
        IGVFeature feature =  codec.decode(line);
        assertEquals("strong GGGCGGGTGGGGCGGG", feature.getName());
    }

    @Test
    public void testSpaceDelimited() throws Exception {
        String line = "chr1 1051161 1051177 strong_GGGCGGGTGGGGCGGG 0.81";
        IGVBEDCodec codec = new IGVBEDCodec();
        IGVFeature feature =  codec.decode(line);
        assertEquals("strong_GGGCGGGTGGGGCGGG", feature.getName());
    }

}
