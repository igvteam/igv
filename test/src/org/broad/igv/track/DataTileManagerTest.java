/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
 * DataTileManagerTest.java
 * JUnit based test
 *
 * Created on November 13, 2007, 10:29 PM
 */

package org.broad.igv.track;


/**
 * @author jrobinso
 */
public class DataTileManagerTest {
    /*
   File h5File = new File("chr6.txt.h5");;
   String trackName = "chr6.txt.h5";
   String chr = "chr6";

   public DataTileManagerTest(String testName) {
       super(testName);
   }

   protected void setUp() throws Exception {

   }

   protected void tearDown() throws Exception {
       H5.H5close();
   }

   public void testConstructor() {

       int zoom = 2;

       HDF5TileManager mgr = new HDF5TileManager(h5File.getAbsolutePath(),
               trackName, chr);

       System.out.println("Done");
   }

   public void testGetRawTile() throws Exception {
       int zoom = 2;
       HDF5TileManager mgr = new HDF5TileManager(h5File.getAbsolutePath(),
               trackName, chr);

       HDF5TileManager.RawTile tile = mgr.getRawTile(0);
       System.out.println(tile.startLocation);
       System.out.println(tile.endLocation);
   }

   public void testGetTile() {


       int zoom = 2;
       HDF5TileManager mgr = new HDF5TileManager(h5File.getAbsolutePath(),
               trackName, chr);

       // Location within tile 0
       float loc1 = 8998242.0f;
       DataTile tile0 = mgr.getTileForLocation(zoom, loc1);
       assertEquals(0, tile0.getTileNumber());
       assertEquals(0.01872588f, tile0.getScore(0), 1.0e-9);
       assertEquals(3230822.8f, tile0.getBaseline(0), 1.0e-9);
       assertEquals(417, tile0.getSize());

       // Location within tile 3
       float loc3 = 1.49498976E8f;
       DataTile tile3 = mgr.getTileForLocation(zoom, loc3);
       assertEquals(3, tile3.getTileNumber());
       assertEquals(1.12170968E8f, tile3.getBaseline(0), 1.0e-9);

       // Get range of tiles
       DataTile [] tiles = mgr.getTilesForRange(zoom, loc1, loc3);
       assertEquals(4, tiles.length);
   }
   * */

}
