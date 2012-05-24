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

package org.broad.igv.variant.vcf;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broadinstitute.sting.utils.codecs.vcf.StandardVCFWriter;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.junit.Ignore;
import org.junit.Test;

import java.io.File;

/**
 * User: jacob
 * Date: 2012/05/11
 */
// TODO We don't have the necessary dependencies
@Ignore
public class VCFWriterTest extends AbstractHeadlessTest {

    @Test
    public void testWriteHeader() throws Exception {
        String inpath = TestUtils.DATA_DIR + "vcf/SRP32_v4.0.vcf";
        File outFile = new File(TestUtils.DATA_DIR, "testwriter.vcf");
        VCFWriter writer = new StandardVCFWriter(outFile, null);

        FeatureCodec codec = CodecFactory.getCodec(inpath, genome);
        AbstractFeatureReader<Feature> bfs = AbstractFeatureReader.getFeatureReader(inpath, codec, false);
        Iterable<Feature> iter = bfs.iterator();
        Object header = bfs.getHeader();

        writer.writeHeader((VCFHeader) header);
        writer.close();
    }

    @Test
    public void testWriteRecords() throws Exception {
        String inpath = TestUtils.DATA_DIR + "vcf/SRP32_v4.0.vcf";
        File outFile = new File(TestUtils.DATA_DIR, "testwriter.vcf");
        VCFWriter writer = new StandardVCFWriter(outFile, null);

        FeatureCodec codec = CodecFactory.getCodec(inpath, genome);
        AbstractFeatureReader<Feature> bfs = AbstractFeatureReader.getFeatureReader(inpath, codec, false);
        Iterable<Feature> iter = bfs.iterator();
        Object header = bfs.getHeader();

        writer.writeHeader((VCFHeader) header);

        for (Feature feat : iter) {
            if (feat instanceof VCFVariant) {
                VCFVariant var = (VCFVariant) feat;
                writer.add(var.getVariantContext());
            }
        }

        writer.close();


    }
}
