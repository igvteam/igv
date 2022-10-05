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

package org.broad.igv.feature.tribble.reader;

import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.FileInputStream;
import java.util.List;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 * Date: 1/23/14
 * Time: 4:24 PM
 */
public class BamIndexTest {

    @Test
    public void testParse() throws Exception {

        String path = TestUtils.DATA_DIR + "tabix/refGene.hg19.bed.gz.tbi";

        FileInputStream fis = new FileInputStream(path);
        byte[] bytes = fis.readAllBytes();

        BamIndex tabixIndex = new BamIndex();
        tabixIndex.parse(bytes, null);
        assertNotNull(tabixIndex.sequenceIndexMap);


    }

    @Test
    public void testTBIIndex() throws Exception {

        String path = TestUtils.DATA_DIR + "tabix/refGene.hg19.bed.gz.tbi";
        int refID = 26;
        int beg = 55369194;
        int end = 55369443;

        FileInputStream fis = new FileInputStream(path);
        byte[] bytes = fis.readAllBytes();

        BamIndex tbiIndex = new BamIndex();
        tbiIndex.parse(bytes, null);

        List<BamIndex.Chunk> chunks = tbiIndex.chunksForRange(refID, beg, end);
        assertNotNull(tbiIndex.sequenceIndexMap);
        assertEquals(chunks.size(), 1);
        assertEquals(1640062, chunks.get(0).left.block);

    }

    @Test
    public void testCSIIndex() throws Exception {

        String path = TestUtils.DATA_DIR + "tabix/csi-test.vcf.gz.csi";
        int refID = 0;
        int beg = 1226000;
        int end = 1227000;

        FileInputStream fis = new FileInputStream(path);
        byte[] bytes = fis.readAllBytes();

        CSIIndex csiIndex = new CSIIndex();
        csiIndex.parse(bytes, null);

        List<CSIIndex.Chunk> chunks = csiIndex.chunksForRange(refID, beg, end);
        assertTrue(chunks.size() > 0);

    }

    @Test
    public void testCSIQuery() throws Exception {

        String path = TestUtils.DATA_DIR + "tabix/csi-test.vcf.gz.csi";
        int refID = 0;
        int beg = 1226000;
        int end = 1227000;

        FileInputStream fis = new FileInputStream(path);
        byte[] bytes = fis.readAllBytes();

        CSIIndex csiIndex = new CSIIndex();
        csiIndex.parse(bytes, null);

        List<CSIIndex.Chunk> chunks = csiIndex.chunksForRange(refID, beg, end);
        assertEquals(chunks.size(),  1);
        assertEquals(74500, chunks.get(0).left.block);
        assertEquals(61246, chunks.get(0).left.offset);
        assertEquals(94804, chunks.get(0).right.block);
    }



//
//    test("CSI query - vcf", async function () {
//
//        const chr = "chr1",
//                beg = 1226000,
//                end = 1227000;
//
//        const reader = new FeatureFileReader({
//                format: "vcf",
//                url: require.resolve("./data/tabix/csi-test.vcf.gz"),
//                indexURL: require.resolve("./data/tabix/csi-test.vcf.gz.csi")
//        });
//        await reader.readHeader();
//        const features = await reader.readFeatures(chr, beg, end);
//        assert.ok(features);
//
//        // Verify features are in query range.  The features immediately before and after query range are returned by design
//        const len = features.length;
//        assert.ok(len > 0);
//        for (let i = 1; i < len - 1; i++) {
//            const f = features[i];
//            assert.ok(f.chr === chr && f.end >= beg && f.start <= end);
//        }
//    })
//
//    test("CSI query - gtf", async function () {
//
//        this.timeout(200000);
//
//        const chr = "10",
//                beg = 400000,
//                end = 500000;
//
//        const reader = new FeatureFileReader({
//                format: "gff3",
//                url: "https://s3.amazonaws.com/igv.org.genomes/hg38/Homo_sapiens.GRCh38.94.chr.gff3.gz",
//                indexURL: "https://s3.amazonaws.com/igv.org.genomes/hg38/Homo_sapiens.GRCh38.94.chr.gff3.gz.csi"
//        });
//        await reader.readHeader();
//        const features = await reader.readFeatures(chr, beg, end);
//        assert.ok(features);
//
//        const len = features.length;
//        assert.ok(len > 0);
//        for (let i = 1; i < len - 1; i++) {
//            const f = features[i];
//            assert.ok(f.chr === chr);
//        }
//    })

}
