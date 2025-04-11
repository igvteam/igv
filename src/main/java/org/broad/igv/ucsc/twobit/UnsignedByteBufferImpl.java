package org.broad.igv.ucsc.twobit;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

public class UnsignedByteBufferImpl implements UnsignedByteBuffer {

    ByteBuffer wrappedBuffer;

    public static UnsignedByteBufferImpl loadBinaryBuffer(String path, ByteOrder byteOrder, long start, int size) throws IOException {
        try (SeekableStream is = IGVSeekableStreamFactory.getInstance().getStreamFor(path)) {
            return getUnsignedByteBuffer(is, byteOrder, start, size);
        }
    }

    public static UnsignedByteBufferImpl getUnsignedByteBuffer(SeekableStream is, ByteOrder byteOrder, long start, int size ) throws IOException {
        ByteBuffer bb = ByteBuffer.allocate(size);
        bb.order(byteOrder);
        byte[] bytes = bb.array();
        is.seek(start);
        is.readFully(bytes);
        return new UnsignedByteBufferImpl(bb);
    }

    public static UnsignedByteBufferImpl wrap(byte[] bytes, ByteOrder byteOrder) {
        ByteBuffer bb = ByteBuffer.wrap(bytes);
        bb.order(byteOrder);
        return new UnsignedByteBufferImpl(bb);
    }

    public static UnsignedByteBuffer wrap(byte[] bytes) {
        return wrap(bytes, ByteOrder.LITTLE_ENDIAN);
    }

    private UnsignedByteBufferImpl(ByteBuffer wrappedBuffer) {
        this.wrappedBuffer = wrappedBuffer;
    }

    @Override
    public byte get() {
        return wrappedBuffer.get();
    }

    @Override
    public short getShort() {
        return wrappedBuffer.getShort();
    }

    @Override
    public int getUShort() {
        return Short.toUnsignedInt(wrappedBuffer.getShort());
    }

    @Override
    public int getInt() {
        return wrappedBuffer.getInt();
    }

    @Override
    public long getUInt() {
        return Integer.toUnsignedLong(wrappedBuffer.getInt());
    }

    @Override
    public long getLong() {
        return wrappedBuffer.getLong();
    }

    @Override
    public float getFloat() {
        return wrappedBuffer.getFloat();
    }

    @Override
    public double getDouble() {
        return wrappedBuffer.getDouble();
    }

    @Override
    public byte[] getBytes(int length) {
        byte[] bytes = new byte[length];
        wrappedBuffer.get(wrappedBuffer.position(), bytes);
        wrappedBuffer.position(wrappedBuffer.position() + length);
        return bytes;
    }


    /**
     * Return a null (0) terminated string
     * @return
     */
    @Override
    public String getString() {
        ByteArrayOutputStream bis = new ByteArrayOutputStream(1000);
        int b;
        while ((b = wrappedBuffer.get()) != 0) {
            bis.write((byte) b);
        }
        return new String(bis.toByteArray());
    }

    @Override
    public String getFixedLengthString(int length) {
        byte[] bytes = new byte[length];
        int nonPaddedLength = 0;
        wrappedBuffer.get(bytes);
        for (int i = 0; i < length; i++) {
            if (bytes[i] == 0) break;
            nonPaddedLength++;
        }
        return new String(bytes, 0, nonPaddedLength);
    }

    @Override
    public int position() {
        return wrappedBuffer.position();
    }

    @Override
    public void position(int i) {
        wrappedBuffer.position(i);
    }

    public byte[] array() {
        return wrappedBuffer.array();
    }

    public int remaining() {
        return wrappedBuffer.remaining();
    }

}
