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
import org.apache.log4j.Logger;
import org.broad.igv.feature.AminoAcidManager;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.action.SearchCommand;
import org.broad.igv.ui.event.ViewChange;
import org.broad.igv.ui.panel.FrameManager;
import org.fest.swing.fixture.FrameFixture;
import org.fest.swing.fixture.JTextComponentFixture;
import org.junit.*;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Test the Command bar UI.
 * <p/>
 * User: jacob
 * Date: 2012/05/14
 */
public class IGVCommandBarTest extends AbstractHeadedTest {

    private static Logger log = Logger.getLogger(IGVCommandBarTest.class);

    private static FrameFixture frame;

    //We use asynchronous events, and there are asserts in there
    //Here we check that the events were actually dispatched and handled
    private int expectedEvents = 0;
    private int actualEvents = expectedEvents;

    private boolean registered;

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
    }

    @Before
    public void setUp() throws Exception{
        super.setUp();
        this.actualEvents = 0;
        this.expectedEvents = 0;
        this.registered = false;
    }

    @After
    public void tearDown() throws Exception{
        System.out.println("Expected Events: " + expectedEvents + " Actual Events: " + actualEvents);
        assertTrue("Event handler not triggered properly", actualEvents >= expectedEvents);
        if(registered){
            FrameManager.getDefaultFrame().getEventBus().unregister(this);
        }
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
    private volatile String enterText = null;

    private void tstChromoNav(String chromoText) throws Exception {
        registerEventHandler(1);

        JTextComponentFixture searchFixture = frame.textBox("searchTextField");
        searchFixture.deleteText();
        this.enterText = chromoText;

        //Make sure search box has focus
        searchFixture.focus();
        searchFixture.requireFocused();
        searchFixture.requireEmpty();

        searchFixture.enterText(chromoText);
        frame.button("goButton").click();
    }

    private void registerEventHandler(int expectedEvents) {
        this.expectedEvents += expectedEvents;
        if(!registered){
            FrameManager.getDefaultFrame().getEventBus().register(this);
            this.registered = true;
        }
    }

    @Subscribe
    public void receiveChromoChange(ViewChange.ChromosomeChangeCause e){
        if(e.source instanceof SearchCommand){
            actualEvents++;
            try {
                assertEquals(this.enterText, e.chrName);
            } catch (Exception e1) {
                log.error(e1.getMessage(), e1);
                actualEvents = Integer.MIN_VALUE;
            }
        }else{
            actualEvents = Integer.MIN_VALUE;
            throw new AssertionError("Got a ChromosomeChangeCause event from unexpected source: " + e.source);
        }
    }

    @Test
    public void testChromoNav_CodonTable() throws Exception {
        //Make sure that translation table DOES NOT change when we change chromosomes
        Assume.assumeTrue("hg18.unittest".equals(GenomeManager.getInstance().getGenomeId()));

        tstChromoNav("chr1");
        assertEquals(1, AminoAcidManager.getInstance().getCodonTable().getId());

        tstChromoNav("chrM");
        assertEquals(1, AminoAcidManager.getInstance().getCodonTable().getId());
    }

}
