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

package org.broad.igv.cli_plugin;

import org.broad.igv.data.DataTile;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Exon;
import org.broad.igv.sam.AlignmentDataManager;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.track.GFFFeatureSource;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.Feature;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;

import static org.junit.Assert.*;

/**
 * User: jacob
 * Date: 2013-Feb-15
 */
public class CufflinksPluginTest extends AbstractPluginTest{

    @BeforeClass
    public static void setUpClass() throws Exception {
        pluginPath = "resources/cufflinks_plugin.xml";
        AbstractPluginTest.setUpClass();
    }

    /**
     * Load some alignment data and run cufflinks on it.
     *
     * @param chr
     * @param start
     * @param end
     * @param inFile
     * @return The transcripts.gtf file cufflinks calculated
     */
    private List<Object> basicRunCufflinks(String chr, int start, int end, File inFile) throws Exception{

        //Load some alignment data
        ResourceLocator locator = new ResourceLocator(inFile.getAbsolutePath());
        AlignmentDataManager dataManager = new AlignmentDataManager(locator, genome);
        AlignmentTrack alignmentTrack = new AlignmentTrack(locator, dataManager, genome);
        dataManager.loadAlignments(chr, start, end, null, null);

        PluginSpecReader.Command command = tool.commandList.get(0);

        List<Argument> argumentList = command.argumentList;

        LinkedHashMap<Argument, Object> arguments = new LinkedHashMap<Argument, Object>(argumentList.size());

        //Fill in the input arguments
        for(int argNum = 0; argNum < argumentList.size(); argNum++) {
            Argument argument = argumentList.get(argNum);
            String argName = argument.getName();
            Object value = argument.getDefaultValue();
            if(argName.equalsIgnoreCase("Output Dir")){
                value = TestUtils.TMP_OUTPUT_DIR;
            }else if(argName.equalsIgnoreCase("Track")){
                value = alignmentTrack;
            }

            if(argument.getType() == Argument.InputType.BOOL){
                value = Boolean.parseBoolean((String) value);
            }

            arguments.put(argument, value);
        }

        List<String> commands = Arrays.asList(toolPath);

        PluginFeatureSource featureSource = new PluginFeatureSource(commands, arguments, command.outputList.get(0), pluginPath);

        PluginDataSource geneTrackingSource = new PluginDataSource(genome, commands, arguments, command.outputList.get(1), pluginPath);

        return Arrays.asList(featureSource, geneTrackingSource);
    }

    @Test
    public void testBasic() throws Exception{

        String chr = "chr1";
        int start = 151666493 - 1;
        int end = start + 10000;
        File inFile = new File(TestUtils.LARGE_DATA_DIR, "HG00171.hg18.bam");
        List sources = basicRunCufflinks(chr, start, end, inFile);
        PluginFeatureSource featureSource = (PluginFeatureSource) sources.get(0);
        PluginDataSource geneTrackingSource = (PluginDataSource) sources.get(1);

        Iterator<Feature> features = featureSource.getFeatures(chr, start, end);
        int count = countNonNullFeatures(features);
        assertTrue("No features read" , count > 0);

        assertTrue("No longest feature", geneTrackingSource.getLongestFeature(chr) > 0);
        DataTile dataTile = geneTrackingSource.getRawData(chr, start, end);
        assertNotNull("No fpkm data returned", dataTile);
        assertFalse(dataTile.isEmpty());

    }

    private int countNonNullFeatures(Iterator<Feature> featureIterator){
        int count = 0;
        while(featureIterator.hasNext()){
            Feature feat = featureIterator.next();
            assertNotNull(feat);
            count++;
        }
        return count;
    }

    @Test
    public void testCufflinksExample_test_data() throws Exception{

        String chr = "test_chromosome";
        int start = 0;
        int end = start + 400;
        File inFile = new File(TestUtils.DATA_DIR + "sam", "cufflinks_test_data.sam");
        TestUtils.createIndex(inFile.getAbsolutePath());

        Iterator<Feature> features = ((PluginFeatureSource) basicRunCufflinks(chr, start, end, inFile).get(0)).getFeatures(chr, start, end);

        List<Feature> combinedFeatures = new GFFFeatureSource.GFFCombiner().addFeatures(features).combineFeatures();

        assertEquals("Wrong number of transcripts", 1, combinedFeatures.size());

        BasicFeature transcript = (BasicFeature) combinedFeatures.get(0);

        assertEquals(53 - 1, transcript.getStart());
        assertEquals(550, transcript.getEnd());


        assertEquals(3, transcript.getExonCount());
        String transFPKM = transcript.getAttributes().get("FPKM");
        assertTrue("Transcript FPKM invalid: " + transFPKM, transFPKM.startsWith("105"));

        int[] exonStarts = new int[]{53,351,501};
        int[] exonEnds = new int[]{250,400,550};

        for(int ii= 0; ii < 3; ii++){
            Exon exon = transcript.getExons().get(ii);
            assertEquals(exonStarts[ii] - 1, exon.getStart());
            assertEquals(exonEnds[ii], exon.getEnd());
            String exonFPKM = transcript.getAttributes().get("FPKM");
            assertTrue("Exon FPKM invalid: " + exonFPKM, exonFPKM.startsWith("105"));
            assertTrue(exonFPKM.startsWith("105"));
        }

    }
}
