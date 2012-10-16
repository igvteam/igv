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

import org.broad.igv.feature.AminoAcidManager;
import org.broad.igv.feature.genome.GenomeManager;
import org.fest.swing.fixture.FrameFixture;
import org.fest.swing.fixture.JComboBoxFixture;
import org.fest.swing.fixture.JTextComponentFixture;
import org.junit.Assume;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Test the Command bar UI.
 * <p/>
 * User: jacob
 * Date: 2012/05/14
 */
public class CommandBarTest extends AbstractHeadedTest {

    private static FrameFixture frame;

    /**
     * Because loading IGV takes so long, we only load it once per class.
     * We reset the session in between tests
     *
     * @throws Exception
     */
    @BeforeClass
    public static void setUpClass() throws Exception {
        AbstractHeadedTest.setUpClass();
        frame = new FrameFixture(IGV.getMainFrame());


        //FrameFixture frame = new FrameFixture(IGV.getMainFrame());
        //JPanelFixture contentFixture = frame.panel("contentPane");
        //JPanelFixture commandBar = frame.panel("igvCommandBar");
        //JComboBoxFixture chromoBox = frame.comboBox("chromosomeComboBox");
    }

    /**
     * Basic test showing usage of FEST and checking combo box
     */
    @Test
    public void testChromoBoxContents() throws Exception {
        String[] chromos = frame.comboBox("chromosomeComboBox").contents();
        assertEquals(26, chromos.length);
    }

    @Test
    public void testChromoNav() throws Exception {
        tstChromoNav("chr1");
        tstChromoNav("chr20");
    }

    private void tstChromoNav(String chromoText) throws Exception {
        JTextComponentFixture searchFixture = frame.textBox("searchTextField");
        searchFixture.deleteText();
        String enterText = chromoText;

        //Make sure search box has focus
        searchFixture.focus();
        searchFixture.requireFocused();
        assertEquals("", searchFixture.text());

        searchFixture.enterText(enterText);
        frame.button("goButton").click();

        JComboBoxFixture comboBox = frame.comboBox("chromosomeComboBox");
        comboBox.requireSelection(enterText);
    }

    @Test
    public void testChromoNav_CodonTable() throws Exception {
        //Make sure that translation table DOES NOT change when we change chromosomes
        Assume.assumeTrue("hg18".equals(GenomeManager.getInstance().getGenomeId()));

        tstChromoNav("chr1");
        assertEquals(1, AminoAcidManager.getInstance().getCodonTable().getId());

        tstChromoNav("chrM");
        assertEquals(1, AminoAcidManager.getInstance().getCodonTable().getId());


    }

}
