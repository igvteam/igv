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

import org.broad.igv.util.StringUtils;

import java.io.IOException;
import java.nio.ByteBuffer;

/**
 * @author jrobinso
 */
public class TDFBedTile implements TDFTile {

    int tileStart;
    int[] start;
    int[] end;
    float[][] data;
    String[] name;

    public TDFBedTile(ByteBuffer byteBuffer, int nSamples, TDFTile.Type type) throws IOException {
        this.fill(byteBuffer, nSamples, type);
    }

    public TDFBedTile(int tileStart, int[] start, int[] end, float[][] data) {
        this.tileStart = tileStart;
        this.start = start;
        this.end = end;
        this.data = data;
    }

    public TDFBedTile(int tileStart, int[] start, int[] end, float[][] data, String[] name) {
        this(tileStart, start, end, data);
        this.name = name;
    }

    public int getSize() {
        return start.length;
    }

    public int getTileStart() {
        return tileStart;
    }

    public int getTileEnd() {
        return getSize() == 0 ? 0 : getEndPosition(getSize() - 1);
    }

    public int getStartPosition(int idx) {
        return start[idx];
    }

    public int getEndPosition(int idx) {
        return end[idx];
    }

    public String getName(int idx) {
        return name == null ? null : name[idx];
    }

    public float getValue(int row, int idx) {
        return data[row][idx];
    }

    public void writeTo(BufferedByteWriter fos) throws IOException {

        // File type
        TDFTile.Type type = name == null ? TDFTile.Type.bed : TDFTile.Type.bedWithName;
        fos.putNullTerminatedString(type.toString());

        int nPositions = start.length;
        int nSamples = data.length;


        fos.putInt(nPositions);

        for (int i = 0; i < nPositions; i++) {
            fos.putInt(start[i]);
        }
        for (int i = 0; i < nPositions; i++) {
            fos.putInt(end[i]);
        }

        fos.putInt(nSamples);
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                fos.putFloat(data[i][j]);
            }
        }

        // Optionally record feature names
        if (type == TDFTile.Type.bedWithName) {
            for (int i = 0; i < nPositions; i++) {
                fos.putNullTerminatedString(name[i]);
            }
        }


    }

    private void fill(ByteBuffer byteBuffer, int nSamples, TDFTile.Type type) throws IOException {

        int nPositions = byteBuffer.getInt();
        start = new int[nPositions];
        for (int i = 0; i < nPositions; i++) {
            start[i] = byteBuffer.getInt();
        }
        end = new int[nPositions];
        for (int i = 0; i < nPositions; i++) {
            end[i] = byteBuffer.getInt();
            ;
        }

        int nS = byteBuffer.getInt();
        ;
        assert (nS == nSamples);

        data = new float[nS][nPositions];
        for (int row = 0; row < nS; row++) {
            data[row] = new float[nPositions];
            for (int i = 0; i < nPositions; i++) {
                data[row][i] = byteBuffer.getFloat();
            }
        }

        // Optionally read feature names
        if (type == TDFTile.Type.bedWithName) {
            name = new String[nPositions];
            for (int i = 0; i < nPositions; i++) {
                name[i] = StringUtils.readString(byteBuffer);
            }

        }

    }
}
