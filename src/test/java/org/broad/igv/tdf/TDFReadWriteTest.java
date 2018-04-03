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
