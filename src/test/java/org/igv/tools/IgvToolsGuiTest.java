package org.igv.tools;

import org.igv.AbstractHeadlessTest;
import org.igv.util.FileUtils;
import org.igv.util.TestUtils;
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
