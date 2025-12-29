/*
 * ByteConverter.java
 *
 * Created on August 15, 2007, 8:48 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package org.broad.igv.util;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.IOException;

/**
 * @author jrobinso
 */
public class ByteConverter {

    public static float[] convertBytesToFloats(byte[] bytes) {
        try {
            int nFloats = bytes.length / (Float.SIZE / 8);
            float[] floats = new float[nFloats];
            DataInputStream dis = new DataInputStream(new ByteArrayInputStream(bytes));
            for (int i = 0; i < nFloats; i++) {
                if (dis.available() >= (Float.SIZE / 8)) {
                    floats[i] = dis.readFloat();
                }
            }
            dis.close();
            return floats;
        } catch (IOException ex) {
            ex.printStackTrace();
            throw new RuntimeException(ex);
        }
    }

    public static int[] convertBytesToInts(byte[] bytes) {
        try {
            int nInts = bytes.length / (Integer.SIZE / 8);
            int[] ints = new int[nInts];
            DataInputStream dis = new DataInputStream(new ByteArrayInputStream(bytes));
            for (int i = 0; i < nInts; i++) {
                if (dis.available() >= (Integer.SIZE / 8)) {
                    ints[i] = dis.readInt();
                }
            }
            dis.close();
            return ints;
        } catch (IOException ex) {
            ex.printStackTrace();
            throw new RuntimeException(ex);
        }
    }

    public static long[] convertBytesToLongs(byte[] bytes) {
        try {
            int nLongs = bytes.length / (Long.SIZE / 8);
            long[] longs = new long[nLongs];
            DataInputStream dis = new DataInputStream(new ByteArrayInputStream(bytes));
            for (int i = 0; i < nLongs; i++) {
                if (dis.available() >= (Long.SIZE / 8)) {
                    longs[i] = dis.readLong();
                }
            }
            dis.close();
            return longs;
        } catch (IOException ex) {
            ex.printStackTrace();
            throw new RuntimeException(ex);
        }
    }

    /**
     * Read an array of strings from the data stream. It is assumed that the file contains fixed length ascii arrays
     * of nChars characters each.
     */
    public static String[] convertBytesToStrings(byte[] bytes, int nChars) {
        try {

            int nStrings = (int) (bytes.length / (nChars * (Character.SIZE / 8)));
            String[] strings = new String[nStrings];
            byte[] stringBytes = new byte[nChars * (Character.SIZE / 8)];

            DataInputStream dis = new DataInputStream(new ByteArrayInputStream(bytes));
            for (int i = 0; i < nStrings; i++) {
                if (dis.available() >= (nChars * (Character.SIZE / 8))) {
                    dis.read(stringBytes);
                    strings[i] = new String(stringBytes, "UTF-16BE").trim();
                }
            }
            return strings;

        } catch (IOException ex) {
            ex.printStackTrace();
            throw new RuntimeException(ex);
        }
    }

}
