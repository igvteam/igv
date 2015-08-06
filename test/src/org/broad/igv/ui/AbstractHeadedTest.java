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

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.Track;
import org.broad.igv.util.TestUtils;
import org.junit.*;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.io.IOException;

import static org.junit.Assert.assertTrue;

/**
 * Class from which headed tests can inherit
 * User: jacob
 * Date: 2012/05/24
 */
@Ignore
public class AbstractHeadedTest {

    protected static IGV igv;
    protected static Genome genome;

    @Rule
    public TestRule testTimeout = new Timeout((int) 60000);

    @BeforeClass
    public static void setUpClass() throws Exception {
        assumeNotHeadless();
        TestUtils.setUpTestEnvironment();
        igv = startGUI();

        TestUtils.setAllNames(igv, true);
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
        TestUtils.clearOutputDir();
        stopGUI();
        igv = null;
    }

    @Before
    public void setUp() throws Exception {
        igv.newSession();
        IGV.getMainFrame().requestFocus();

        TestUtils.resetPrefsFile();
        TestUtils.resetTestUserDefinedGenomes();
        IGV.getInstance().getContentPane().getCommandBar().refreshGenomeListComboBox();
    }

    @After
    public void tearDown() throws Exception {
        TestUtils.resetPrefsFile();
        TestUtils.resetTestUserDefinedGenomes();
    }


    /**
     * Start GUI with default genome file
     *
     * @return
     * @throws java.io.IOException
     */
    public static IGV startGUI() throws IOException {
        return startGUI(TestUtils.defaultGenome);
    }


    /**
     * Load a gui with the specified genome file.
     * No genome is loaded if null
     *
     * @param genomeFile
     * @return
     * @throws IOException
     */
    protected static IGV startGUI(String genomeFile) throws IOException {
        Globals.setHeadless(false);
        IGV igv;
        //If IGV is already open, we get the instance.
        if (IGV.hasInstance()) {
            igv = IGV.getInstance();
            IGV.getMainFrame().setVisible(true);
            System.out.println("Using old IGV");
        } else {
            JFrame frame = new JFrame();
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            Main.open(frame);
            System.out.println("Started new IGV");
            igv = IGV.getInstance();
            assertTrue(IGV.getInstance().waitForNotify(1000));
        }
        if (genomeFile != null) {
            igv.loadGenome(genomeFile, null, true);
            genome = igv.getGenomeManager().getCurrentGenome();
        }
        return igv;
    }


    /**
     * This closes the IGV window.
     */
    public static void stopGUI() {
        if (!IGV.hasInstance()) {
            return;
        }

        IGV.getMainFrame().setVisible(false);
        IGV.getMainFrame().dispose();
        IGV.destroyInstance();
    }

    public static void assumeNotHeadless() {
        boolean headless = true;
        try {
            GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
            headless = ge.isHeadless();
        } catch (Exception e) {
            e.printStackTrace();
        } catch (Error e) {
            //Really not sure why this ever happens, maybe just jenkins issues
            e.printStackTrace();
        }
        if (headless) {
            System.out.println("You are trying to start a GUI in a headless environment. Aborting test");
        }
        org.junit.Assume.assumeTrue(!headless);
    }


    public void testTest() throws Exception {
        java.util.List<Track> tracks = IGV.getInstance().getAllTracks();
        System.out.println("# tracks: " + tracks.size());
        for (Track track : tracks) {
            System.out.println(track.getName());
        }

        java.util.List<Track> featureTracks = IGV.getInstance().getTrackPanel(IGV.FEATURE_PANEL_NAME).getTracks();
        System.out.println(featureTracks.size());
    }

    protected static String rewriteRestoreSession(String sessionPath) throws Exception{
        sessionPath = (TestUtils.replaceTestPaths(new File(sessionPath))).getAbsolutePath();
        IGV.getInstance().doRestoreSession(sessionPath, null, false);
        return sessionPath;
    }
}
