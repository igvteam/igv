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
        String fname = "Unigene.sample.bed";
        File baseInputFile = new File(TestUtils.DATA_DIR + "bed", fname);
        File inputFile = new File(TestUtils.DATA_DIR + "folder.with.periods", fname);
        inputFile.delete();
        inputFile.deleteOnExit();
        FileUtils.copyFile(baseInputFile, inputFile);


        String output = IgvToolsGui.getDefaultOutputText(inputFile.getAbsolutePath(), IgvToolsGui.Tool.SORT);
        File outputFile = new File(output);
        outputFile.delete();
        outputFile.deleteOnExit();

        //The file doesn't actually exist, but we should
        //be able to write to the path
        assertFalse(outputFile.exists());
        FileWriter writer = new FileWriter(outputFile);
        assertNotNull(writer);
    }
}
