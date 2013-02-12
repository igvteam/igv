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

import com.google.common.eventbus.Subscribe;
import org.broad.igv.feature.AminoAcidManager;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.action.SearchCommand;
import org.broad.igv.ui.event.ViewChange;
import org.broad.igv.ui.panel.FrameManager;
import org.fest.swing.fixture.FrameFixture;
import org.fest.swing.fixture.JComboBoxFixture;
import org.fest.swing.fixture.JTextComponentFixture;
import org.junit.*;

import javax.swing.*;

import static org.junit.Assert.assertEquals;

/**
 * Test the Command bar UI.
 * <p/>
 * User: jacob
 * Date: 2012/05/14
 */
public class IGVCommandBarTest extends AbstractHeadedTest {

    private static FrameFixture frame;

    //We use asynchronous events, and there are asserts in there
    //Here we check that the events were actually dispatched and handled
    private int expectedEvents = 0;
    private int actualEvents = expectedEvents;

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

    @Before
    public void setUp() throws Exception{
        super.setUp();
        this.actualEvents = 0;
        this.expectedEvents = 0;
    }

    @After
    public void tearDown() throws Exception{
        assertEquals("Event handler not triggered properly", expectedEvents, actualEvents);
        super.tearDown();
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

    //Tacky state variable, but eh, it's just a test
    private String enterText = null;

    private void tstChromoNav(String chromoText) throws Exception {
        registerEventHandler(2);

        JTextComponentFixture searchFixture = frame.textBox("searchTextField");
        searchFixture.deleteText();
        enterText = chromoText;

        //Make sure search box has focus
        searchFixture.focus();
        searchFixture.requireFocused();
        assertEquals("", searchFixture.text());

        searchFixture.enterText(enterText);
        frame.button("goButton").click();
    }

    private void registerEventHandler(int expectedEvents) {
        this.expectedEvents += expectedEvents;
        FrameManager.getDefaultFrame().getEventBus().register(this);
    }

    /**
     * This is a little funny, we are actually responding to a spurious ChromosomeChangeCause
     * sent when the dropdown changes
     * @param e
     */
    @Subscribe
    public void receiveChromoChange(ViewChange.ChromosomeChangeCause e){
        if(e.source instanceof JComboBox){
            actualEvents++;
            JComboBoxFixture comboBox = frame.comboBox("chromosomeComboBox");
            comboBox.requireSelection(enterText);
        }else if(e.source instanceof SearchCommand){
            actualEvents++;
            assertEquals(enterText, e.chrName);
        }else{
            actualEvents = Integer.MIN_VALUE;
            throw new AssertionError("Got a ChromosomeChangeCause event from unexpected source: " + e.source);
        }
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
