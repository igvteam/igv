/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.tdf;

import java.io.ByteArrayOutputStream;
import java.io.IOException;

/*
 * 
 */
public class BufferedByteWriter {

    /**
     * Reports the total number of bytes written
     */

    ByteArrayOutputStream buffer;


    public BufferedByteWriter() {
        this(8192);
    }


    public BufferedByteWriter(int size) {
        if (size <= 0) {
            throw new IllegalArgumentException("Buffer size <= 0");
        }
        buffer = new ByteArrayOutputStream(size);
    }

    public byte [] getBytes() {
        return buffer.toByteArray();
    }

    public int bytesWritten() {
        return buffer.size();
    }

    /**
     * Writes the specified byte to this buffered output stream.
     *
     * @param b the byte to be written.
     * @throws IOException if an I/O error occurs.
     */
    public void put(int b) throws IOException {
        buffer.write(b);
    }

    public void put(byte b[]) throws IOException {
        buffer.write(b);
    }

    /**
     * Writes <code>len</code> bytes from the specified byte array
     * starting at offset <code>off</code> to this buffered output stream.
     *
     * @param b   the data.
     * @param off the start offset in the data.
     * @param len the number of bytes to write.
     * @throws IOException if an I/O error occurs.
     */
    public void put(byte b[], int off, int len) throws IOException {
        buffer.write(b, off, len);
    }

    public void putInt(int v) throws IOException {
        buffer.write((v >>> 0) & 0xFF);
        buffer.write((v >>> 8) & 0xFF);
        buffer.write((v >>> 16) & 0xFF);
        buffer.write((v >>> 24) & 0xFF);
    }


    public void putShort(short v) throws IOException {
        buffer.write((v >>> 0) & 0xFF);
        buffer.write((v >>> 8) & 0xFF);
    }


    public void putFloat(float f) throws IOException {
        int v = Float.floatToIntBits(f);
        putInt(v);
    }

    public void putDouble(Double f) throws IOException {
        long v = Double.doubleToLongBits(f);
        putLong(v);
    }


    /**
     * Writes a <code>long</code> to the underlying output stream as eight
     * bytes, little endian. In no exception is thrown, the counter
     * <code>written</code> is incremented by <code>8</code>.
     *
     * @param v a <code>long</code> to be written.
     * @throws IOException if an I/O error occurs.
     * @see java.io.FilterOutputStream#out
     */
    public void putLong(long v) throws IOException {
        buffer.write((byte) (v >>> 0));
        buffer.write((byte) (v >>> 8));
        buffer.write((byte) (v >>> 16));
        buffer.write((byte) (v >>> 24));
        buffer.write((byte) (v >>> 32));
        buffer.write((byte) (v >>> 40));
        buffer.write((byte) (v >>> 48));
        buffer.write((byte) (v >>> 56));
    }


    public void putNullTerminatedString(String s) throws IOException {
        buffer.write(s.getBytes());
        buffer.write((byte) 0);
    }

}
