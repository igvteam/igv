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

package org.broad.igv.data;

import org.apache.commons.lang.ArrayUtils;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.util.TestUtils;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureCodec;
import org.junit.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static junit.framework.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2013-Mar-11
 */
public class GenomeSummaryDataTest extends AbstractHeadlessTest {

    private static final String testDataPath = TestUtils.DATA_DIR + "bed/chr1-2_peaks.bed";

    @Override
    public void setUp() throws Exception {
        super.setUp();
        TestUtils.createIndex(testDataPath);
    }

    private Iterable getFeatures(String chr) throws Exception {
        FeatureCodec codec = CodecFactory.getCodec(testDataPath, genome);
        AbstractFeatureReader bfs = AbstractFeatureReader.getFeatureReader(testDataPath, codec, true);
        if(chr == null) return bfs.iterator();
        return bfs.query(chr, 0, Integer.MAX_VALUE);
    }

    /**
     * Build GenomeSummaryData from sorted (by chromosome and start position) set of features
     * @param genome
     * @param trackId
     * @param features
     * @return
     */
    private GenomeSummaryData buildGenomeSummaryData(Genome genome, String trackId, Iterable<BasicFeature> features, double scale){
        GenomeSummaryData genomeSummaryData = new GenomeSummaryData(genome, new String[]{trackId});
        if(scale > 0){
            genomeSummaryData.setScale(scale);
        }

        List<Integer> startLocations = new ArrayList<Integer>();
        List<Float> data = new ArrayList<Float>();
        String lastChr = null;

        for(BasicFeature feature: features){
            String chr = feature.getChr();
            //Finish off last chromosome
            if(lastChr != null && !chr.equals(lastChr)){
                Map<String, float[]> dMap = new HashMap<String, float[]>();
                dMap.put(trackId, ArrayUtils.toPrimitive(data.toArray(new Float[data.size()])));
                genomeSummaryData.addData(lastChr, ArrayUtils.toPrimitive(startLocations.toArray(new Integer[startLocations.size()])), dMap);
                startLocations.clear();
                data.clear();
            }

            startLocations.add(feature.getStart());
            data.add(feature.getScore());

            lastChr = chr;
        }

        Map<String, float[]> dMap = new HashMap<String, float[]>();
        dMap.put(trackId, ArrayUtils.toPrimitive(data.toArray(new Float[data.size()])));
        genomeSummaryData.addData(lastChr, ArrayUtils.toPrimitive(startLocations.toArray(new Integer[startLocations.size()])), dMap);

        return genomeSummaryData;
    }

    @Test
    public void testSinglePointsSingleChromo_01() throws Exception{
        tstAccumulate("chr1", 6, 100.0);
    }

    @Test
    public void testSinglePointsSingleChromo_02() throws Exception{
        tstAccumulate("chr2", 8, 10.0);
    }

    @Test
    public void testCombinePointsSingleChromo_01() throws Exception{
        tstAccumulate("chr1", 3, 1000.0);
    }

    @Test
    public void testCombinePointsSingleChromo_02() throws Exception{
        tstAccumulate("chr2", 5, 1000.0);
    }

    @Test
    public void testSinglePointsMultiChromo() throws Exception{
        tstAccumulate(null, 14, 10.0);
    }

    @Test
    public void testCombinePointsMultiChromo_01() throws Exception{
        tstAccumulate(null, 13, 100.0);
    }

    @Test
    public void testCombinePointsMultiChromo_02() throws Exception{
        tstAccumulate(null, 7, -1);
    }


    private void tstAccumulate(String chr, int expLocations, double scale) throws Exception{
        GenomeSummaryData genomeSummaryData = buildGenomeSummaryData(genome, testDataPath, getFeatures(chr), scale);

        assertEquals(0, genomeSummaryData.skippedChromosomes.size());

        assertEquals(expLocations, genomeSummaryData.nDataPts);
        assertEquals(expLocations, genomeSummaryData.getLocations().length);
        assertEquals(expLocations, genomeSummaryData.getData(testDataPath).length);
    }
}
