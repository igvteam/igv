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
