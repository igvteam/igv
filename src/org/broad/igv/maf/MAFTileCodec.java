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
package org.broad.igv.maf;

import org.broad.igv.maf.MAFTile.MASequence;

import java.io.DataOutputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 */
public class MAFTileCodec {

    /**
     * Serialize an MATile to a stream
     */
    public void encode(MAFTile maTile, DataOutputStream os) {

        try {
            os.writeInt(maTile.start);
            os.writeInt(maTile.end);
            os.writeInt(maTile.alignedSequences.size());
            for (Map.Entry<String, MASequence> entry : maTile.alignedSequences.entrySet()) {
                os.writeUTF(entry.getKey());
                os.writeUTF(entry.getValue().bases);
            }

            os.writeInt(maTile.gapAdjustedIdx.length);
            for (int i = 0; i < maTile.gapAdjustedIdx.length; i++) {
                os.writeInt(maTile.gapAdjustedIdx[i]);
            }
            os.flush();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public MAFTile decode(java.io.DataInputStream inputStream) {

        try {

            int start = inputStream.readInt();
            int end = inputStream.readInt();
            int sz = inputStream.readInt();
            if (sz == 0) {
                return new MAFTile(start, end);
            }
            Map<String, String> bases = new HashMap(sz);
            for (int i = 0; i < sz; i++) {
                String speciesId = inputStream.readUTF();
                String b = inputStream.readUTF();
                bases.put(speciesId, b);
                //   tile.alignedSequences.put(speciesId, new MASequence(tile, bases));
            }

            sz = inputStream.readInt();
            int[] gapAdjustedCoordinates = new int[sz];
            for (int i = 0; i < sz; i++) {
                gapAdjustedCoordinates[i] = inputStream.readInt();
            }
            return new MAFTile(start, end, bases, gapAdjustedCoordinates);

        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
    }

}
