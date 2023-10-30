package org.broad.igv.util;

import htsjdk.samtools.seekablestream.SeekableStream;

import java.io.ByteArrayOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

public class UnsignedByteBuffer {

    ByteBuffer wrappedBuffer;

    public static UnsignedByteBuffer loadBinaryBuffer(SeekableStream is, ByteOrder byteOrder, long start, int size) throws IOException {
        ByteBuffer bb = ByteBuffer.allocate(size);
        bb.order(byteOrder);
        byte[] bytes = bb.array();
        is.seek(start);
        is.readFully(bytes);
        return new UnsignedByteBuffer(bb);
    }


    public UnsignedByteBuffer(ByteBuffer wrappedBuffer) {
        this.wrappedBuffer = wrappedBuffer;
    }

    public byte get() {
        return wrappedBuffer.get();
    }

    public short getShort() {
        return wrappedBuffer.getShort();
    }

    public int getUShort() {
        return Short.toUnsignedInt(wrappedBuffer.getShort());
    }

    public int getInt() {
        return wrappedBuffer.getInt();
    }

    public long getUInt() {
        return Integer.toUnsignedLong(wrappedBuffer.getInt());
    }

    public long getLong() {
        return wrappedBuffer.getLong();
    }

    public String getString() throws IOException {
        ByteArrayOutputStream bis = new ByteArrayOutputStream(1000);
        int b;
        while ((b = wrappedBuffer.get()) != 0) {
            if(b < 0) {
                throw new EOFException();
            }
            bis.write((byte) b);
        }
        return new String(bis.toByteArray());
    }

    public String getFixedLengthString(int length) throws IOException {

        byte [] bytes = new byte[length];
        for(int i=0; i<length; i++) {
            bytes[i] = wrappedBuffer.get();
        }
        return new String(bytes);

    }

    public void position(int i) {
        wrappedBuffer.position(i);
    }

    public byte[] array() {
        return wrappedBuffer.array();
    }
}
