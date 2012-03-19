/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.tdf;

import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.assertEquals;

/**
 * File file, TrackType trackType, String trackLine, String[] trackNames
 *
 * @author jrobinso
 */
public class TDFReadWriteTest {


    static String[] trackNames = {"sample 1", "sample 2", "sample 3"};
    static TrackType type = TrackType.COPY_NUMBER;
    static String trackLine = "track line";
    static List<WindowFunction> wfs = Arrays.asList();

    public TDFReadWriteTest() {


    }

    @BeforeClass
    public static void setUpClass() throws Exception {

    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }


    /**
     * Test of createGroup method, of class IBFWriter.
     */
    @Test
    public void testCreateGroup() {
        File testFile = new File("test1.tdf");
        testFile.deleteOnExit();

        String groupName = "/chr1";

        TDFWriter writer = new TDFWriter(testFile, "hg18", type, trackLine, trackNames, wfs, false);
        TDFGroup group = writer.createGroup(groupName);
        group.setAttribute("attKey", "attValue");
        group.setAttribute("attKey2", "attValue2");
        writer.closeFile();

        TDFReader reader = TDFReader.getReader(testFile.getAbsolutePath());
        group = reader.getGroup(groupName);
        assertEquals(groupName, group.getName());
        assertEquals(2, group.attributes.size());
        assertEquals("attValue", group.attributes.get("attKey"));
        assertEquals("attValue2", group.attributes.get("attKey2"));
        reader.close();
    }

    /**
     * Test of createDataset method, of class IBFWriter.
     */
    @Test
    public void testCreateFile() {
        File testFile = new File("test2.tdf");
        testFile.deleteOnExit();

        TDFWriter writer = new TDFWriter(testFile, "hg18", type, trackLine, trackNames, wfs, false);
        writer.closeFile();

        TDFReader reader = TDFReader.getReader(testFile.getAbsolutePath());
        String[] tn = reader.getTrackNames();
        assertEquals(trackNames.length, tn.length);
        for (int i = 0; i < tn.length; i++) {
            assertEquals(trackNames[i], tn[i]);
        }
        assertEquals(type, reader.getTrackType());
        assertEquals(trackLine, reader.getTrackLine());
        reader.close();


    }

    /**
     * Test of writeTile method, of class IBFWriter.
     */
    @Test
    public void testWriteTile() throws IOException {
        writeTile("test3.tdf", false);
    }


    /**
     * Test of writeTile method, of class IBFWriter.
     */
    @Test
    public void testWriteCompressedTile() throws IOException {
        writeTile("test4.tdf", true);
    }


    public void writeTile(String file, boolean gzipped) throws IOException {

        File testFile = new File(file);
        testFile.deleteOnExit();

        String dsName = "/chr1/zoom1/medians";
        int tileNumber = 0;
        int start = 10;
        float span = 1;


        float[][] data = new float[trackNames.length][100000];
        for (int i = 0; i < trackNames.length; i++) {
            for (int j = 0; j < 10; j++) {
                data[i][j] = (float) Math.random();
            }
        }
        //System.out.println("----------------");

        TDFWriter writer = new TDFWriter(testFile, "hg18", type, trackLine, trackNames, wfs, gzipped);
        TDFDataset dataset = writer.createDataset(dsName, TDFDataset.DataType.FLOAT, 1000, 1);
        dataset.setAttribute("attr1", "value1");
        dataset.setAttribute("attr2", "value2");

        TDFFixedTile tile = new TDFFixedTile(start, start, span, data);
        writer.writeTile(dsName, tileNumber, tile);
        writer.closeFile();

        TDFReader reader = TDFReader.getReader(testFile.getAbsolutePath());

        TDFDataset ds = reader.getDataset(dsName);

        assertEquals("value1", ds.getAttribute("attr1"));
        assertEquals("value2", ds.getAttribute("attr2"));

        TDFTile dsTile = reader.readTile(ds, 0);

        for (int i = 0; i < trackNames.length; i++) {
            assertEquals(start + i * span, dsTile.getStartPosition(i), 1.0e-6);
            for (int j = 0; j < 10; j++) {
                assertEquals(data[i][j], dsTile.getValue(i, j), 1.0e-6);
            }
        }
        reader.close();
    }


    public static void main(String[] args) {
        org.junit.runner.JUnitCore.runClasses(TDFReadWriteTest.class);

    }
}
