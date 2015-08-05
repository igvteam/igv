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
