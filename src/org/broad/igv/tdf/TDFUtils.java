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

import org.broad.igv.Globals;
import org.broad.igv.track.WindowFunction;

import java.io.BufferedOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Map;

/**
 * @author jrobinso
 */
public class TDFUtils {

    public static void main(String[] args) throws FileNotFoundException {
        //chr1:170,118,422-170,764,140
        TDFUtils.dumpIndex("/Users/jrobinso/IGV/gbm/GBM.transcriptome__agilentg4502a_07_2__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.data.tdf");
        //TDFUtils.dumpAllTiles("/Users/jrobinso/projects/client/compressed.tdf");
        //TDFUtils.tdfToBedgraph("/Users/jrobinso/IGV/time_course/cebpb_0.merged.bam.tdf",
        //        "/Users/jrobinso/IGV/time_course/test.wig");
        //TDFUtils.dumpTile("/Users/jrobinso/IGV/allaml.dataset.gct.tdf", "/chr4/raw", 130);
        //TDFUtils.dumpDatasets("/Users/jrobinso/IGV/affy_1552717.gct.tdf");
        //TDFUtils.dumpDatasets("/Volumes/igv/tools/OV-1116.capture.tumor.20b.cov.tdf");
    }

    public static void dumpIndex(String tdfFile) {

        TDFReader reader = TDFReader.getReader(tdfFile);
        reader.dumpIndex();
        reader.close();;
    }

    public static void tdfToBedgraph(String tdfFile, String bedGraphFile) throws FileNotFoundException {

        TDFReader reader = null;
        PrintStream ps = null;

        try {
            reader = TDFReader.getReader(tdfFile);
            ps = new PrintStream(new BufferedOutputStream(new FileOutputStream(bedGraphFile)));

            String trackLine = reader.getTrackLine();
            if (trackLine != null && trackLine.length() > 0) {
                ps.println(trackLine);
            }
            for (String dsName : reader.getDatasetNames()) {

                String[] tokens = dsName.split("/");
                String chrName = tokens[1];
                if (!chrName.equals(Globals.CHR_ALL) && dsName.contains("raw")) {

                    TDFDataset ds = reader.getDataset(dsName);
                    for (int i = 0; i < ds.nTiles; i++) {
                        TDFTile tile = ds.getTile(i);
                        if (tile != null) {
                            dumpTileData(reader, chrName, tile, ps);
                        }
                    }
                }
            }
        } finally {
            if (reader != null) reader.close();
            if (ps != null) ps.close();
        }


    }


    public static void dumpRootAttributes(String ibfFile) {
        TDFReader reader = TDFReader.getReader(ibfFile);
        System.out.println("Track line = " + reader.getTrackLine());
        TDFGroup group = reader.getGroup("/");
        for (Map.Entry<String, String> entries : group.attributes.entrySet()) {
            System.out.println(entries.getKey() + " = " + entries.getValue());
        }

        System.out.println(reader.getTrackLine());
    }

    public static void dumpDatasets(String ibfFile) {
        TDFReader reader = TDFReader.getReader(ibfFile);
        System.out.println("DATASETS");
        for (String dsName : reader.getDatasetNames()) {
            System.out.println(dsName);
            TDFDataset ds = reader.getDataset(dsName);

            System.out.println("Attributes");
            for (Map.Entry<String, String> entry : ds.attributes.entrySet()) {
                System.out.println("\t" + entry.getKey() + " = " + entry.getValue());
            }
            System.out.println();

            System.out.println("Tile Positions");
            for (int i = 0; i < ds.nTiles; i++) {
                System.out.print("\t" + ds.tilePositions[i]);
            }
            System.out.println();

        }
    }

    public static void dumpAllTiles(String ibfFile) {
        TDFReader reader = TDFReader.getReader(ibfFile);
        System.out.println("DATASETS");
        for (String dsName : reader.getDatasetNames()) {
            System.out.println(dsName);
            TDFDataset ds = reader.getDataset(dsName);

            for (int i = 0; i < ds.nTiles; i++) {
                TDFTile tile = ds.getTile(i);
                if (tile != null) {
                    System.out.println("Tile: " + i);
                    dumpTileData(reader, "", tile, System.out);
                }
            }
        }
    }

    public static void dumpTile(String ibfFile, String dsName, int tileNumber) {

        TDFReader reader = TDFReader.getReader(ibfFile);
        TDFDataset ds = reader.getDataset(dsName);
        TDFTile tile = reader.readTile(ds, tileNumber);
        if (tile == null) {
            System.out.println("Null tile: " + dsName + " [" + tileNumber + "]");
        } else {
            dumpTileData(reader, "", tile, System.out);
        }
    }

    private static void dumpTileData(TDFReader reader, String chrName, TDFTile tile, PrintStream ps) {
        int nTracks = reader.getTrackNames().length;
        int nBins = tile.getSize();
        if (nBins > 0) {
            for (int b = 0; b < nBins; b++) {
                ps.print(chrName);
                ps.print("\t" + tile.getStartPosition(b));
                ps.print("\t" + tile.getEndPosition(b));
                for (int t = 0; t < nTracks; t++) {
                    ps.print("\t" + tile.getValue(t, b));
                }
                ps.println();
            }
        }
    }

    public static void dumpRange(String ibfFile, String dsName, int startLocation, int endLocation) {
        TDFReader reader = TDFReader.getReader(ibfFile);
        TDFDataset ds = reader.getDataset(dsName);

        int tileWidth = ds.tileWidth;
        int startTile = startLocation / tileWidth;
        int endTile = endLocation / tileWidth;

        for (int tileNumber = startTile; tileNumber <= endTile; tileNumber++) {
            TDFTile tile = reader.readTile(ds, tileNumber);
            if (tile == null) {
                System.out.println("Null tile: " + dsName + " [" + tileNumber + "]");
            } else {
                int nTracks = reader.getTrackNames().length;
                int nBins = tile.getSize();
                if (nBins > 0) {
                    for (int b = 0; b < nBins; b++) {
                        int start = tile.getStartPosition(b);
                        int end = tile.getEndPosition(b);
                        if (start > endLocation) {
                            break;
                        }
                        if (end >= startLocation) {
                            System.out.print(tile.getStartPosition(b));
                            for (int t = 0; t < nTracks; t++) {
                                System.out.print("\t" + tile.getValue(t, b));
                            }
                            System.out.println();
                        }
                    }
                }
            }

        }
    }


    /*
     *     magic number (4 bytes)
    version
    index position
    index size (bytes)
    # of window functions
    [window functions]
    track type (string)
    track line (string)
    # of tracks
    [track names]
     */

    public static void dumpSummary(String ibfFile) {

        TDFReader reader = TDFReader.getReader(ibfFile);

        System.out.println("Version: " + reader.getVersion());
        System.out.println("Window Functions");
        for (WindowFunction wf : reader.getWindowFunctions()) {
            System.out.println("\t" + wf.toString());
        }

        System.out.println("Tracks");
        String[] trackNames = reader.getTrackNames();
        for (String trackName : trackNames) {
            System.out.println(trackName);
        }
        System.out.println();

        System.out.println("DATASETS");
        for (String dsName : reader.getDatasetNames()) {
            System.out.println(dsName);
            TDFDataset ds = reader.getDataset(dsName);

            System.out.println("Attributes");
            for (Map.Entry<String, String> entry : ds.attributes.entrySet()) {
                System.out.println("\t" + entry.getKey() + " = " + entry.getValue());
            }
            System.out.println();

            System.out.println("Tiles");

            int nTracks = trackNames.length;
            int tracksToShow = Math.min(4, nTracks);

            for (int i = 0; i < ds.nTiles; i++) {
                TDFTile tile = reader.readTile(ds, i);
                if (tile != null) {
                    System.out.print("  " + i);
                    /*int nBins = tile.getSize();
                    int binsToShow = Math.min(4, nBins);
                    for (int b = 0; b < binsToShow; b++) {
                        System.out.print(tile.getStartPosition(b));
                        for (int t = 0; t < tracksToShow; t++) {
                            float value = tile.getValue(0, b);
                            if (!Float.isNaN(value)) {
                                System.out.print("\t" + tile.getValue(t, b));
                            }
                        }
                        System.out.println();
                    }
                    */
                }
            }
            System.out.println();
            System.out.println();
        }

        System.out.println("GROUPS");
        for (String name : reader.getGroupNames()) {
            System.out.println(name);
            TDFGroup group = reader.getGroup(name);

            System.out.println("Attributes");
            for (Map.Entry<String, String> entry : group.attributes.entrySet()) {
                System.out.println("\t" + entry.getKey() + " = " + entry.getValue());
            }
            System.out.println();

        }
    }

    /*
    chr1:241356465-241356657
    chr1:241358198-241358223
    chr1:241359291-241359329
    
    chr4:119691730-119691768
    chr4:119692843-119692868
    chr4:119694419-119694611

     */
}
