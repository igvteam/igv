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

package org.broad.igv.data;

import org.apache.commons.lang.ArrayUtils;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.FeatureCodec;
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
