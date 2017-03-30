package org.broad.igv.sam.lite;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

/**
 * Created by jrobinso on 3/10/17.
 */
public class BGUnzip {
    public static final int BGZIP_HEADER_LENGTH = 18;

    // Uncompress data,  assumed to be series of bgzipped blocks

    public static byte[] blockUnzip(byte[] data) throws IOException {

        Inflater inflater = new Inflater(true);

        int ptr = 0;

        int lim = data.length - BGZIP_HEADER_LENGTH;
        ByteArrayOutputStream outputStream = new ByteArrayOutputStream(data.length);

        while (ptr < lim) {

            int xlen =  unpackInt16(data, ptr + 10);
            int si1 = data[ptr + 12];
            int si2 = data[ptr + 13];
            int slen =  unpackInt16(data, 14);
            int bsize = unpackInt16(data, ptr + 16) + 1;

            int start = BGZIP_HEADER_LENGTH + ptr;    // Start of CDATA
            int remainder = data.length - start;
            if (remainder < (bsize + 8)) break;

            int uncLength = unpackInt32(data, ptr + bsize - 4);

            inflater.reset();
            inflater.setInput(data, start, remainder);

            byte [] output = new byte[uncLength];
            final int inflatedBytes;
            try {
                inflatedBytes = inflater.inflate(output, 0, uncLength);
            } catch (DataFormatException e) {
                // Can happen near end of file;
                break;
            }

            outputStream.write(output, 0, inflatedBytes);

            ptr += bsize;    // Advance to next block

        }

        return outputStream.toByteArray();
    }

    private static int unpackInt16(final byte[] buffer, final int offset) {
        return ((buffer[offset] & 0xFF) |
                ((buffer[offset+1] & 0xFF) << 8));
    }


    private static int unpackInt32(final byte[] buffer, final int offset) {
        return ((buffer[offset] & 0xFF) |
                ((buffer[offset + 1] & 0xFF) << 8) |
                ((buffer[offset + 2] & 0xFF) << 16) |
                ((buffer[offset + 3] & 0xFF) << 24));
    }
}
