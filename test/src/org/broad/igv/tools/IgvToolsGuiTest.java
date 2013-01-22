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

import org.broad.igv.ui.AbstractHeadedTest;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;
import java.io.FileWriter;

import static junit.framework.Assert.assertFalse;
import static junit.framework.Assert.assertNotNull;

/**
 * Test of the IGVToolsGui. Technically don't need to extend
 * AbstractHeadedTest, since we don't need an IGV instance, but
 * it doesn't hurt
 * User: jacob
 * Date: 2013-Jan-22
 */
public class IgvToolsGuiTest extends AbstractHeadedTest{

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

        //Create IGVTools window, initialize params for this test
        IgvToolsGui mainWindow = new IgvToolsGui();
        mainWindow.setModal(false);
        mainWindow.setVisible(true);
        mainWindow.setTool(IgvToolsGui.SORT);

        mainWindow.setInputFieldText(inputFile.getAbsolutePath());
        String output = mainWindow.getOutputFieldText();
        File outputFile = new File(output);

        //The file doesn't actually exist, but we should
        //be able to write to the path
        System.out.println(outputFile);
        assertFalse(outputFile.exists());
        FileWriter writer = new FileWriter(outputFile);
        assertNotNull(writer);
    }
}
