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
