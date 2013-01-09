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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.tdf;

import org.broad.igv.util.StringUtils;
import org.broad.igv.util.collections.LRUCache;

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

    DataType dataType;
    int tileWidth;
    long[] tilePositions;  // File position in TDF file
    int[] tileSizes;       // Tile size in bytes
    int nTiles;
    LRUCache<String, TDFTile> cache = new LRUCache(this, 20);
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

    public List<TDFTile> getTiles() {
        List<TDFTile> tiles = new ArrayList<TDFTile>();
        for (int t = 0; t < nTiles; t++) {
            TDFTile tile = getTile(t);
            if (tile != null) {
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

    public void clearCache() {
        cache.clear();
    }


    public DataType getDataType() {
        return dataType;
    }

    public int getTileWidth() {
        return tileWidth;
    }


}
