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

package org.broad.igv.ui;

import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.StringUtils;
import org.broad.igv.util.TestUtils;
import org.junit.*;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import java.io.File;
import java.util.Collection;

import static org.junit.Assert.*;

/**
 * Test of class main. In general will use this to see that IGV
 * starts properly, given startup parameters. Since we will be starting
 * IGV in a non-standard way, do NOT inherit from AbstractHeadedTest
 * User: jacob
 * Date: 2012/05/04
 */
public class MainTest {

    @Rule
    public TestRule testTimeout = new Timeout((int) 1e5);

    @BeforeClass
    public static void setUpClass() throws Exception {
        AbstractHeadedTest.assumeNotHeadless();
        TestUtils.setUpTestEnvironment();
    }

    @Before
    public void setUp() throws Exception {
        AbstractHeadedTest.stopGUI();
    }


    @After
    public void tearDown() throws Exception {
        AbstractHeadedTest.stopGUI();
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
        AbstractHeadedTest.tearDownClass();
    }


    /**
     * Test that loading IGV with a startup file and
     * locus loads that file and locus
     *
     * @throws Exception
     */
    @Test
    public void testFileLocusArgs() throws Exception {

        String trackName = "NA12878.SLX.sample.bam";
        String filePath = TestUtils.DATA_DIR + "bam/" + trackName;
        String genome = "hg18";
        String locus = "chr1:9,703,210-9,727,225";

        String[] args = new String[]{filePath, locus, "-g", genome};

        //Need to wait for IGV to load file, genome, and move to locus
        IGV igv = startMain(args, 60000);

        //System.out.println(IGV.getInstance());

        assertEquals(genome, igv.getGenomeManager().getGenomeId());
        TestUtils.assertTrackLoaded(igv, trackName);

        String actLocus = FrameManager.getDefaultFrame().getFormattedLocusString();
        assertEquals(locus, actLocus);
    }

    @Test
    public void testFileWithSpaces() throws Exception {
        String trackName = "test.wig";
        String filePath = TestUtils.DATA_DIR + "folder with spaces/" + trackName;
        String[] args = new String[]{filePath};

        //Need to wait for IGV to start and load file
        IGV igv = startMain(args, 10000);

        TestUtils.assertTrackLoaded(igv, trackName);
    }

    @Test
    public void testFileURLWithSpaces() throws Exception {
        String trackName = "test.wig";
        String dir = StringUtils.encodeURL("folder with spaces");
        String absFilePath = (new File(TestUtils.DATA_DIR)).getAbsolutePath();
        absFilePath = absFilePath.replace("\\", "/");
        String filePath = "file://" + absFilePath + "/" + dir + "/" + trackName;
        String[] args = new String[]{filePath};

        //Need to wait for IGV to start and load file
        IGV igv = startMain(args, 10000);

        TestUtils.assertTrackLoaded(igv, trackName);
    }

    @Test
    public void testRemoteURLWithSpaces() throws Exception {
        String trackName = "test.wig";
        String dir = StringUtils.encodeURL("folder with spaces");
        String absFilePath = "www.broadinstitute.org/igvdata/unit_test_files";
        String filePath = String.format("http://%s/%s/%s", absFilePath, dir, trackName);
        String[] args = new String[]{filePath};

        //Need to wait for IGV to start and load file
        IGV igv = startMain(args, 10000);

        TestUtils.assertTrackLoaded(igv, trackName);
    }

    /**
     * Test loading a genome not in the display list , by id
     *
     * @throws Exception
     */
    @Test
    public void testLoadGenomeById() throws Exception {
        String genomeId = "mm7";
        Collection<GenomeListItem> genomeListItems = GenomeManager.getInstance().getGenomes();
        for (GenomeListItem gen : genomeListItems) {
            assertNotSame("Bad test setup, test genome in display list", gen.getId(), genomeId);
        }

        String[] args = new String[]{"-g", genomeId};
        IGV igv = startMain(args, 10000);

        assertEquals(igv.getGenomeManager().getGenomeId(), genomeId);
    }

    /**
     * Test loading a genome not in the display list, by full path
     *
     * @throws Exception
     */
    @Test
    public void testLoadGenomeByPath() throws Exception {
        String genomePath = TestUtils.DATA_DIR + "genomes/canFam2.unittest.genome";
        String genomeId = "canFam2.unittest";

        String[] args = new String[]{"-g", genomePath};
        IGV igv = startMain(args, 10000);

        assertEquals(igv.getGenomeManager().getGenomeId(), genomeId);
    }

    /**
     * Test loading a fasta not in the display list, by full path
     *
     * @throws Exception
     */
    @Test
    public void testLoadFastaByPath() throws Exception {
        String genomePath = TestUtils.DATA_DIR + "fasta/ecoli_out.padded.fasta";
        String genomeId = genomePath;

        String[] args = new String[]{"-g", genomePath};
        IGV igv = startMain(args, 5000);

        assertEquals(igv.getGenomeManager().getGenomeId(), genomeId);
    }

    /**
     * Test loading a fasta in working directory
     * This was GH #66 and IGV-2042
     *
     * @throws Exception
     */
    @Test
    public void testLoadFastaWorkDir() throws Exception {

        String genomeFiName = "ecoli_out.padded.fasta";
        String srcGenomePath = TestUtils.DATA_DIR + "fasta/" + genomeFiName;
        File destFile = new File(genomeFiName);
        destFile.delete();
        destFile.deleteOnExit();
        FileUtils.copyFile(new File(srcGenomePath), destFile);

        String[] args = new String[]{"-g", genomeFiName};
        IGV igv = startMain(args, 5000);

        assertEquals(igv.getGenomeManager().getGenomeId(), genomeFiName);
    }


    /*
     * Example  version string:  2.3.27
     */
    @Test
    public void testVersionComparisons() throws Exception {

        String v1 = "2.3.27";

        Main.Version version1 = Main.Version.getVersion(v1);
        assertEquals(2, version1.getMajor());
        assertEquals(3, version1.getMinor());
        assertEquals(27, version1.getBuild());

        String v2 = "2.3.27";
        Main.Version version2 = Main.Version.getVersion(v2);
        assertEquals(2, version2.getMajor());
        assertEquals(3, version2.getMinor());
        assertEquals(27, version2.getBuild());

        assertFalse(version1.lessThan(version2));
        assertFalse(version2.lessThan(version1));

        String v3 = "2.3.26";
        Main.Version version3 = Main.Version.getVersion(v3);
        assertFalse(version1.lessThan(version3));

        String v4 = "2.3.28";
        Main.Version version4 = Main.Version.getVersion(v4);
        assertTrue(version1.lessThan(version4));

    }

    private IGV startMain(String[] args, long timeout) {
        Main.main(args);
        IGV igv = IGV.getInstance();
        assertTrue(igv.waitForNotify(timeout));
        return igv;
    }


}
