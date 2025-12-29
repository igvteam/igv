package org.broad.igv.util;

import java.io.ByteArrayOutputStream;
import java.io.EOFException;
import java.io.FilterInputStream;
import java.io.IOException;
import java.io.InputStream;

public class LittleEndianInputStream extends htsjdk.tribble.util.LittleEndianInputStream {
    byte[] buffer = new byte[8];

    public LittleEndianInputStream(InputStream in) {
        super(in);
    }

    public int readUShort() throws IOException {

        int byte1 = in.read();
        int byte2 = in.read();
        if(byte2 < 0) {
            throw new EOFException();
        } else {
            return ((byte2 << 24 >>> 16) + (byte1 << 24 >>> 24));
        }
    }

    public void readFully(byte[] b) throws IOException {
        readFully(b, 0, b.length);
    }

    public void readFully(byte[] b, int off, int len) throws IOException {
        if(len < 0) {
            throw new IndexOutOfBoundsException("len is negative");
        } else {
            int total;
            int result;
            for(total = 0; total < len; total += result) {
                result = in.read(b, off + total, len - total);
                if(result == -1) {
                    break;
                }
            }

        }
    }

    @Override
    public String readString() throws IOException {
        ByteArrayOutputStream bis = new ByteArrayOutputStream(1000);
        int b;
        while ((b = in.read()) != 0) {
            if(b < 0) {
                throw new EOFException();
            }
            bis.write((byte) b);
        }
        return new String(bis.toByteArray());
    }
}
