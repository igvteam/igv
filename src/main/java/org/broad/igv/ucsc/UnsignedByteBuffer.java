package org.broad.igv.ucsc;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

public class UnsignedByteBuffer {

    ByteBuffer wrappedBuffer;


    public static UnsignedByteBuffer loadBinaryBuffer(String path, ByteOrder byteOrder, long start, int size) throws IOException {
        try(SeekableStream is = IGVSeekableStreamFactory.getInstance().getStreamFor(path)) {
            ByteBuffer bb = ByteBuffer.allocate(size);
            bb.order(byteOrder);
            byte[] bytes = bb.array();
            is.seek(start);
            is.readFully(bytes);
            return new UnsignedByteBuffer(bb);
        }
    }

    public static UnsignedByteBuffer wrap(byte[] bytes, ByteOrder byteOrder) {
        ByteBuffer bb = ByteBuffer.wrap(bytes);
        bb.order(byteOrder);
        return new UnsignedByteBuffer(bb);
    }


    private UnsignedByteBuffer(ByteBuffer wrappedBuffer) {
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

    public double getDouble() {
        return wrappedBuffer.getDouble();
    }

    public String getString()  {
        ByteArrayOutputStream bis = new ByteArrayOutputStream(1000);
        int b;
        while ((b = wrappedBuffer.get()) != 0) {
            if(b < 0) {
                return new String(bis.toByteArray());
            }
            bis.write((byte) b);
        }
        return new String(bis.toByteArray());
    }

    public String getFixedLengthString(int length)  {

        byte [] bytes = new byte[length];
        for(int i=0; i<length; i++) {
            bytes[i] = wrappedBuffer.get();
        }
        return new String(bytes);

    }
    public long position() {
        return wrappedBuffer.position();
    }
    public void position(int i) {
        wrappedBuffer.position(i);
    }

    public byte[] array() {
        return wrappedBuffer.array();
    }

    public long remaining() {
        return wrappedBuffer.remaining();
    }


}
