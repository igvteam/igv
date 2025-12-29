package org.broad.igv.tools.parsers;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.tools.converters.GCTtoIGVConverter;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

/**
 * @author jrobinso
 * @date Oct 9, 2010
 */
public class GCTtoIGVConverterTest extends AbstractHeadlessTest {

    @Test
    public void testDescriptionMapping() throws IOException {
        String gctFile = TestUtils.DATA_DIR + "gct/igv_test2.gct";
        File igvFile = new File(TestUtils.DATA_DIR + "gct/igv_test2.gct.igv");
        igvFile.deleteOnExit();

        ResourceLocator locator = new ResourceLocator(gctFile);
        GCTtoIGVConverter.convert(locator, igvFile, null, 50000, null, genome);
    }


    @Test
    public void testAffyMapping() throws IOException {
        String gctFile = TestUtils.DATA_DIR + "gct/affy_human.gct";
        File igvFile = new File(TestUtils.DATA_DIR + "gct/affy_human.gct.igv");
        igvFile.deleteOnExit();

        ResourceLocator locator = new ResourceLocator(gctFile);
        GCTtoIGVConverter.convert(locator, igvFile, null, 50000, null, genome);
    }

}
