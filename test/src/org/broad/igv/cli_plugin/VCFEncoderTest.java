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

package org.broad.igv.cli_plugin;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.igv.variant.VariantTrack;
import org.broad.igv.variant.vcf.VCFVariant;
import org.broad.igv.variant.vcf.VCFWriterTest;
import org.broad.tribble.Feature;
import org.broadinstitute.variant.vcf.VCFHeader;
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
