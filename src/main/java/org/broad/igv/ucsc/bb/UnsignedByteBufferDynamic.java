package org.broad.igv.ucsc.bb;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.ucsc.twobit.UnsignedByteBuffer;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.ByteArrayOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * A UnsignedByteBuffer that refills its backing buffer as needed from the underlying file resource. This class was
 * created specifically to load the chromTree of a BB file, where the start position is known but total size is
 * not.  It might have general utility but has not been tested with any other use case.
 */
public class UnsignedByteBufferDynamic implements UnsignedByteBuffer {

    /**
     * The current wrapped buffer.  This is updated from the file as needed.
     */
    private ByteBuffer wrappedBuffer;

    /**
     * The file offset corresponding to the start of the current wrapped buffer
     */
    private long offset;

    /**
     * The file offset at object construction.
     */
    private long originalOffset;

    int bufferSize;
    ByteOrder byteOrder;
    String path;

    public static UnsignedByteBufferDynamic loadBinaryBuffer(String path, ByteOrder byteOrder, long offset, int size) throws IOException {
        UnsignedByteBufferDynamic b = new UnsignedByteBufferDynamic(path, byteOrder, offset, size);
        b.updateBuffer();
        return b;
    }

    private UnsignedByteBufferDynamic(String path, ByteOrder byteOrder, long offset, int bufferSize) {
        this.path = path;
        this.byteOrder = byteOrder;
        this.offset = offset;
        this.originalOffset = offset;
        this.bufferSize = bufferSize;
    }

    private void advanceBuffer() {
        offset += (this.wrappedBuffer == null ? 0 : this.wrappedBuffer.position());
        updateBuffer();
    }

    private void updateBuffer() {
        try (SeekableStream is = IGVSeekableStreamFactory.getInstance().getStreamFor(path)) {
            this.wrappedBuffer = ByteBuffer.allocate(bufferSize);
            this.wrappedBuffer.order(byteOrder);
            byte[] bytes = this.wrappedBuffer.array();
            is.seek(offset);
            try {
                is.readFully(bytes);
            } catch (EOFException e) {
                // This can happen if near the end of the file and is o.k., in fact expected
            }
        } catch (IOException e) {
            // TODO ??
        }
    }

    @Override
    public byte get() {
        if (wrappedBuffer.remaining() < 1) {
            advanceBuffer();
        }
        return wrappedBuffer.get();
    }

    @Override
    public short getShort() {
        if (wrappedBuffer.remaining() < 2) {
            advanceBuffer();
        }
        return wrappedBuffer.getShort();
    }

    @Override
    public int getUShort() {
        if (wrappedBuffer.remaining() < 2) {
            advanceBuffer();
        }
        return Short.toUnsignedInt(wrappedBuffer.getShort());
    }

    @Override
    public int getInt() {
        if (wrappedBuffer.remaining() < 4) {
            advanceBuffer();
        }
        return wrappedBuffer.getInt();
    }

    @Override
    public long getUInt() {
        if (wrappedBuffer.remaining() < 4) {
            advanceBuffer();
        }
        return Integer.toUnsignedLong(wrappedBuffer.getInt());
    }

    @Override
    public long getLong() {
        if (wrappedBuffer.remaining() < 8) {
            advanceBuffer();
        }
        return wrappedBuffer.getLong();
    }

    @Override
    public float getFloat() {
        if (wrappedBuffer.remaining() < 4) {
            advanceBuffer();
        }
        return wrappedBuffer.getFloat();
    }

    @Override
    public double getDouble() {
        if (wrappedBuffer.remaining() < 8) {
            advanceBuffer();
        }
        return wrappedBuffer.getDouble();
    }

    @Override
    public byte[] getBytes(int length) {
        if (wrappedBuffer.remaining() < length) {
            advanceBuffer();
        }
        byte[] bytes = new byte[length];
        wrappedBuffer.get(wrappedBuffer.position(), bytes);
        wrappedBuffer.position(wrappedBuffer.position() + length);
        return bytes;
    }


    /**
     * Return a null (0) terminated string.  This method assumes short strings, and will fail if string length is > 1000
     *
     * @return
     */
    @Override
    public String getString() {
        if (wrappedBuffer.remaining() < 1000) {
            advanceBuffer();
        }
        ByteArrayOutputStream bis = new ByteArrayOutputStream(1000);
        int b;
        while ((b = wrappedBuffer.get()) != 0) {
            bis.write((byte) b);
        }
        return new String(bis.toByteArray());
    }

    @Override
    public String getFixedLengthString(int length) {
        if (wrappedBuffer.remaining() < length) {
            advanceBuffer();
        }
        byte[] bytes = new byte[length];
        int nonPaddedLength = 0;
        wrappedBuffer.get(bytes);
        for (int i = 0; i < length; i++) {
            if (bytes[i] == 0) break;
            nonPaddedLength++;
        }
        return new String(bytes, 0, nonPaddedLength);
    }

    /**
     * Position is interpreted as relative to the original file offset
     *
     * @return
     */
    @Override
    public int position() {
        return (int) (offset - originalOffset) + wrappedBuffer.position();
    }

    /**
     * Position is interpreted as relative to the original file offset
     *
     * @return
     */
    @Override
    public void position(int i) {

        /**
         * Position is interpreted relative to the original file position for this buffer
         * @return
         */
        int newBufferPosition = i - (int) (offset - originalOffset);

        if (newBufferPosition < 0 || newBufferPosition > wrappedBuffer.position() + wrappedBuffer.remaining()) {
            offset = originalOffset + i;
            updateBuffer();
        } else {
            wrappedBuffer.position(newBufferPosition);
        }
    }

    @Override
    public int remaining() {
        return wrappedBuffer.remaining();
    }



}
