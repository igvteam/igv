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

package org.broad.igv.cli_plugin;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.igv.variant.VariantTrack;
import org.broad.igv.variant.vcf.VCFVariant;
import org.broad.igv.variant.vcf.VCFWriterTest;
import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFHeader;
import org.junit.Test;

import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static junit.framework.Assert.assertEquals;

/**
 * @author jacob
 * @date 2013-Jul-02
 */
public class VCFEncoderTest extends AbstractHeadlessTest {


    @Test
    public void basicIOTest() throws Exception {
        String inputPath = TestUtils.DATA_DIR + "vcf/hc_mod.vcf";
        TestUtils.createIndex(inputPath);

        VariantTrack inputTrack = (VariantTrack) (new TrackLoader()).load(new ResourceLocator(inputPath), genome).get(0);

        VCFEncoder encoder = new VCFEncoder();
        Argument arg = new Argument(null, Argument.InputType.VARIANT_TRACK, null, null, null, null, false, null);

        Map<Argument, Object> arguments = new HashMap(1);
        arguments.put(arg, inputTrack);
        encoder.setInputs(null, arguments, arg);

        String chr = "chr20";
        int start = 10499005;
        int end = 10499500;

        String outpath = TestUtils.TMP_OUTPUT_DIR + "hc_mod_out.vcf";
        OutputStream outputStream = new FileOutputStream(outpath);

        List<VCFVariant> infeatures = new ArrayList<VCFVariant>();
        for (Feature feat : inputTrack.getFeatures(chr, start, end)) {
            infeatures.add((VCFVariant) feat);
        }

        encoder.encodeAll(outputStream, infeatures.iterator());

        TestUtils.createIndex(outpath);
        VariantTrack outputTrack = (VariantTrack) (new TrackLoader()).load(new ResourceLocator(inputPath), genome).get(0);

        VCFWriterTest.assertHeadersEquals((VCFHeader) inputTrack.getHeader(), (VCFHeader) outputTrack.getHeader());

        int n = 0;
        for (Feature var1 : inputTrack.getFeatures(chr, start, end)) {
            VCFVariant var0 = infeatures.get(n++);
            VCFWriterTest.assertVCFVariantsEqual(var0, (VCFVariant) var1);
        }

        assertEquals("Different number of features from input", infeatures.size(), n);

    }


}
