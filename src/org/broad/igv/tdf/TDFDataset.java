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
    LRUCache<String, TDFTile> cache = new LRUCache(20);
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
