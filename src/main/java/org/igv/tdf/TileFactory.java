package org.igv.tdf;


import org.igv.util.StringUtils;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * @author jrobinso
 */
public class TileFactory {

    public static TDFTile createTile(byte[] buffer, int nSamples) throws IOException {

        ByteBuffer byteBuffer = ByteBuffer.wrap(buffer);
        byteBuffer.order(ByteOrder.LITTLE_ENDIAN);

        String typeString = StringUtils.readString(byteBuffer);
        TDFTile.Type type = TDFTile.Type.valueOf(typeString);

        switch (type) {
            case fixedStep:
                return new TDFFixedTile(byteBuffer, nSamples);
            case variableStep:
                return new TDFVaryTile(byteBuffer, nSamples);
            case bed:
            case bedWithName:
                return new TDFBedTile(byteBuffer, nSamples, type);
            default:
                throw new RuntimeException("Unknown tile type: " + type.toString());
        }
    }


}
