package org.broad.igv.ucsc.twobit;

public interface UnsignedByteBuffer {
    byte get();

    short getShort();

    int getUShort();

    int getInt();

    long getUInt();

    long getLong();

    float getFloat();

    double getDouble();

    String getString();

    String getFixedLengthString(int length);

    int position();

    void position(int i);


    int remaining();

    byte[] getBytes(int valSize);
}
