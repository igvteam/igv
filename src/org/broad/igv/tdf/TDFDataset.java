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
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.tdf;

import org.broad.igv.util.LRUCache;
import org.broad.igv.util.StringUtils;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Represents the data for a particular chromosome and zoom level
 *
 * @author jrobinso
 */
public class TDFDataset extends TDFEntity {

    public enum DataType {

        BYTE, SHORT, INT, FLOAT, DOUBLE, STRING
    }

    ;
    DataType dataType;
    int tileWidth;
    long[] tilePositions;
    int[] tileSizes;
    int nTiles;
    LRUCache<String, TDFTile> cache = new LRUCache(this,20);
    // TODO -- refactor this dependency out
    TDFReader reader;

    public TDFDataset(String name, DataType dataType, int tileWidth, int nTiles) {
        super(name);
        this.dataType = dataType;
        this.tileWidth = tileWidth;
        this.nTiles = nTiles;
        this.tilePositions = new long[nTiles];
        this.tileSizes = new int[nTiles];

        // Initialize tile positions to -1.  This indicates a blank tile and is the default
        Arrays.fill(tilePositions, -1);

    }

    public TDFDataset(String name, ByteBuffer byteBuffer, TDFReader reader) throws IOException {
        super(name);
        this.reader = reader;
        fill(byteBuffer);
    }

    public void write(BufferedByteWriter dos) throws IOException {

        writeAttributes(dos);

        writeString(dos, dataType.toString());
        dos.putFloat(tileWidth);
        //       dos.writeFloat(binWidth);
        dos.putInt(tilePositions.length);
        for (int i = 0; i < tilePositions.length; i++) {
            dos.putLong(tilePositions[i]);
            dos.putInt(tileSizes[i]);
        }

    }

    private void fill(ByteBuffer byteBuffer) throws IOException {

        // Attributes
        readAttributes(byteBuffer);

        String typeString = StringUtils.readString(byteBuffer);
        dataType = TDFDataset.DataType.valueOf(typeString);

        // TODO -- change tileWidth to int ?
        tileWidth = (int) byteBuffer.getFloat();


        nTiles = byteBuffer.getInt();
        tilePositions = new long[nTiles];
        tileSizes = new int[nTiles];

        for (int i = 0; i < nTiles; i++) {
            tilePositions[i] = byteBuffer.getLong();
            tileSizes[i] = byteBuffer.getInt();
        }

    }

    // TODO -- this uses an implied linear index.  Abstract index or replace
    // with general interval index
    // TODO -- gather all non-cached tiles and read in one chunk.  See BAM
    // alignment reader class
    public List<TDFTile> getTiles(int startLocation, int endLocation) {

        List<TDFTile> tiles = new ArrayList();
        int startTile = (int) (startLocation / tileWidth);
        int endTile = (int) (endLocation / tileWidth);
        for (int t = startTile; t <= endTile; t++) {
            TDFTile tile = getTile(t);
            if (tile != null && tile.getSize() > 0) {
                tiles.add(tile);
            }
        }
        return tiles;

    }

    // TDFTile computeTile(TDFDataset ds, int t, List<LocusScore> scores, String chr)
    synchronized TDFTile getTile(int t) {
        String key = getName() + "_" + t;

        TDFTile tile = null;
        if (!cache.containsKey(key)) {
            tile = reader.readTile(this, t);
            cache.put(key, tile);
        } else {
            tile = cache.get(key);
        }
        return tile;
    }

}
