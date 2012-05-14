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

import org.broad.igv.util.TestUtils;
import org.fest.swing.fixture.FrameFixture;
import org.fest.swing.fixture.JComboBoxFixture;
import org.fest.swing.fixture.JPanelFixture;
import org.fest.swing.fixture.JTextComponentFixture;
import org.junit.*;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThat;

/**
 * Test the Command bar UI.
 *
 * User: jacob
 * Date: 2012/05/14
 */
public class CommandBarTest {

    private static IGV igv;
    private static FrameFixture frame;

    /**
     * Because loading IGV takes so long, we only load it once per class.
     * We reset the session in between tests
     * @throws Exception
     */
    @BeforeClass
    public static void setUpClass() throws Exception{
        igv = TestUtils.startGUI();
        frame = new FrameFixture(IGV.getMainFrame());


        //FrameFixture frame = new FrameFixture(IGV.getMainFrame());
        //JPanelFixture contentFixture = frame.panel("contentPane");
        //JPanelFixture commandBar = frame.panel("igvCommandBar");
        //JComboBoxFixture chromoBox = frame.comboBox("chromosomeComboBox");
    }

    @AfterClass
    public static void tearDownClass() throws Exception{
        TestUtils.stopGUI();
        igv = null;
    }

    @Before
    public void setUp() throws Exception {
        igv.resetSession(null);

    }

    @After
    public void tearDown() throws Exception {
        //TODO
    }

    /**
     * Basic test showing usage of FEST and checking combo box
     */
    @Test
    public void testChromoBoxContents() throws Exception{
        String[] chromos = frame.comboBox("chromosomeComboBox").contents();
        assertEquals(26, chromos.length);
    }

    @Test
    public void testChromoNav() throws Exception{
        JTextComponentFixture searchFixture = frame.textBox("searchTextField");
        String enterText = "chr1";
        searchFixture.enterText(enterText);
        frame.button("goButton").click();

        JComboBoxFixture comboBox = frame.comboBox("chromosomeComboBox");
        comboBox.requireSelection(enterText);
    }

}
