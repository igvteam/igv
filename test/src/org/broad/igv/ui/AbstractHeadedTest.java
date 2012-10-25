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

package org.broad.igv.ui;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.TestUtils;
import org.junit.*;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import javax.swing.*;
import java.awt.*;
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
        igv.resetSession(null);
        IGV.getMainFrame().requestFocus();
    }

    @After
    public void tearDown() throws Exception {

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
            igv.loadGenome(genomeFile, null);
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
        }catch(Error e){
            //Really not sure why this ever happens, maybe just jenkins issues
            e.printStackTrace();
        }
        if (headless) {
            System.out.println("You are trying to start a GUI in a headless environment. Aborting test");
        }
        org.junit.Assume.assumeTrue(!headless);
    }
}
