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
public class TDFFixedTile implements TDFTile {

    int tileStart;
    double span;
    int start;
    float[][] data;

    public TDFFixedTile(ByteBuffer byteBuffer, int nSamples) throws IOException {
        this.fill(byteBuffer, nSamples);
    }

    public TDFFixedTile(int tileStart, int start, double span, float[][] data) {
        this.tileStart = tileStart;
        this.span = span;
        this.data = data;
        this.start = start;
    }

    public int getTileStart() {
        return start;
    }

    public int getTileEnd() {
        return getSize() == 0 ? 0 : getEndPosition(getSize() - 1);
    }

    public int getStartPosition(int idx) {
        return start + (int) (idx * span);
    }

    public int getEndPosition(int idx) {
        return start + (int) ((idx + 1) * span);
    }

    public String getName(int idx) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public float getValue(int row, int idx) {
        return data[row][idx];
    }

    public int getSize() {
        return (data == null ? 0 : data[0].length);
    }


    // TODO -- record "type",  extent (longest feature), other stuff

    public void writeTo(BufferedByteWriter fos) throws IOException {

        fos.putNullTerminatedString(TDFTile.Type.fixedStep.toString());
        fos.putInt(getSize());
        fos.putInt(start);
        fos.putFloat((float) span);
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                fos.putFloat(data[i][j]);
            }
        }

    }

    public void fill(ByteBuffer byteBuffer, int nSamples) throws IOException {

        int nPositions = byteBuffer.getInt();
        start = byteBuffer.getInt();
        span = byteBuffer.getFloat();

        data = new float[nSamples][nPositions];
        for (int sample = 0; sample < nSamples; sample++) {
            data[sample] = new float[nPositions];
            for (int i = 0; i < nPositions; i++) {
                data[sample][i] = byteBuffer.getFloat();
            }
        }

    }


    /**
     * This should never be called, but is provided to satisfy the interface
     *
     * @return
     */
    public int[] getStart() {
        int nPts = data[0].length;
        int[] startArray = new int[nPts];
        for (int i = 0; i < nPts; i++) {
            startArray[i] = start + (int) (i * span);
        }
        return startArray;
    }

    public int[] getEnd() {
        int nPts = data[0].length;
        int[] endArray = new int[nPts];
        for (int i = 0; i < nPts; i++) {
            endArray[i] = start + (int) ((i + 1) * span);
        }
        return endArray;
    }

    public float[] getData(int trackNumber) {
        return data[trackNumber];  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String[] getNames() {
        return null;
    }

}
