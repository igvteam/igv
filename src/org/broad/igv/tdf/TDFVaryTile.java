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

import java.io.IOException;
import java.nio.ByteBuffer;

/**
 * @author jrobinso
 */
public class TDFVaryTile implements TDFTile {

    int tileStart;
    double span;
    int[] start;
    float[][] data;

    public TDFVaryTile(ByteBuffer byteBuffer, int nSamples) throws IOException {
        this.fill(byteBuffer, nSamples);
    }

    public TDFVaryTile(int tileStart, double span, int[] start, float[][] data) {
        this.tileStart = tileStart;
        this.span = span;
        this.start = start;
        this.data = data;
    }


    public int getSize() {
        return start.length;
    }

    public int getTileStart() {
        return tileStart;
    }

    public int getStartPosition(int idx) {
        return start[idx];
    }

    public int getEndPosition(int idx) {
        return (int) (start[idx] + span);
    }

    public String getName(int idx) {
        return null;
    }

    public float getValue(int row, int idx) {
        return data[row][idx];
    }

    public void writeTo(BufferedByteWriter fos) throws IOException {

        // File type
        fos.putNullTerminatedString(TDFTile.Type.variableStep.toString());

        fos.putInt(tileStart);
        fos.putFloat((float) span);

        int nPositions = start.length;
        int nSamples = data.length;

        fos.putInt(nPositions);

        for (int i = 0; i < nPositions; i++) {
            fos.putInt(start[i]);
        }

        fos.putInt(nSamples);
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                fos.putFloat(data[i][j]);
            }
        }
    }

    public void fill(ByteBuffer byteBuffer, int nSamples) throws IOException {

        tileStart = byteBuffer.getInt();
        span = byteBuffer.getFloat();

        int nPositions = byteBuffer.getInt();
        start = new int[nPositions];
        for (int i = 0; i < nPositions; i++) {
            start[i] = byteBuffer.getInt();
        }

        int nS = byteBuffer.getInt();
        assert (nS == nSamples);

        data = new float[nS][nPositions];
        for (int row = 0; row < nS; row++) {
            data[row] = new float[nPositions];
            for (int i = 0; i < nPositions; i++) {
                data[row][i] = byteBuffer.getFloat();
            }
        }

    }


    public int[] getStart() {
        return start;
    }

    public int[] getEnd() {
        int [] end = new int[start.length];
        for(int i=0; i<end.length; i++) {
            end[i] = (int) (start[i] + span);
        }
        return end;
    }

    public float[] getData(int trackNumber) {
        return data[trackNumber];  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String[] getNames() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }


}
