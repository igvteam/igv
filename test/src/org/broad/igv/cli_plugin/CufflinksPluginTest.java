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

import org.broad.igv.sam.AlignmentDataManager;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.Feature;
import org.junit.BeforeClass;
import org.junit.Ignore;

import java.io.File;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * User: jacob
 * Date: 2013-Feb-15
 */
@Ignore
public class CufflinksPluginTest extends AbstractPluginTest{

    @BeforeClass
    public static void setUpClass() throws Exception {
        pluginPath = "resources/cufflinks_plugin.xml";
        AbstractPluginTest.setUpClass();
    }

    //@Test
    public void testBasic() throws Exception{


        PluginSpecReader.Command command = tool.commandList.get(0);

        List<Argument> argumentList = command.argumentList;

        //Set output dir
        LinkedHashMap<Argument, Object> arguments = new LinkedHashMap<Argument, Object>(argumentList.size());
        int argnum = 0;
        arguments.put(argumentList.get(argnum++), TestUtils.TMP_OUTPUT_DIR);
        arguments.put(argumentList.get(argnum), argumentList.get(argnum).getDefaultValue());
        argnum++;

        String chr = "chr1";
        int start = 151666493 - 1;
        int end = start + 10000;

        //Load some alignment data
        File inBAMFile = new File(TestUtils.LARGE_DATA_DIR, "HG00171.hg18.bam");
        ResourceLocator locator = new ResourceLocator(inBAMFile.getAbsolutePath());

        AlignmentDataManager dataManager = new AlignmentDataManager(locator, genome);
        AlignmentTrack alignmentTrack = new AlignmentTrack(locator, dataManager, genome);

        dataManager.loadAlignments(chr, start, end, null, null);

        arguments.put(argumentList.get(argnum++), alignmentTrack);

        List<String> commands = Arrays.asList(toolPath);

        PluginFeatureSource combinedFeatureSource = new PluginFeatureSource(commands, arguments, command.parser, pluginPath);

        Iterator<Feature> features = combinedFeatureSource.getFeatures("chr1", start, end);
        while(features.hasNext()){
            System.out.println(features.next().getStart());
        }

    }
}
