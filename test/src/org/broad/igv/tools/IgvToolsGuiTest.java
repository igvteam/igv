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

package org.broad.igv.tools;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;
import java.io.FileWriter;

import static junit.framework.Assert.assertFalse;
import static junit.framework.Assert.assertNotNull;

/**
 * Test of the IGVToolsGui. As of this writing the tests
 * are actually headless.
 * User: jacob
 * Date: 2013-Jan-22
 */
public class IgvToolsGuiTest extends AbstractHeadlessTest{

    /**
     * Copy a file into a folder with periods in the name,
     * check that the default naming functions properly (rather than
     * splitting on the first period we want to split on the last)
     * @throws Exception
     */
    @Test
    public void testSortName_FolderPeriods() throws Exception{
        String fname = "unigene.sample.bed";
        File baseInputFile = new File(TestUtils.DATA_DIR + "bed", fname);
        File inputFile = new File(TestUtils.DATA_DIR + "folder.with.periods", fname);
        inputFile.delete();
        inputFile.deleteOnExit();
        FileUtils.copyFile(baseInputFile, inputFile);


        String output = IgvToolsGui.getDefaultOutputText(inputFile.getAbsolutePath(), IgvToolsGui.Tool.SORT);
        File outputFile = new File(output);
        outputFile.deleteOnExit();

        //The file doesn't actually exist, but we should
        //be able to write to the path
        assertFalse(outputFile.exists());
        FileWriter writer = new FileWriter(outputFile);
        assertNotNull(writer);
    }
}
