/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * MiscTests.java
 *
 * Created on November 16, 2007, 2:48 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.util;

//~--- JDK imports ------------------------------------------------------------

import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author jrobinso
 */
public class MiscTests {

    /**
     * Creates a new instance of MiscTests
     */
    public static void main(String[] args) throws FileNotFoundException, IOException {
        String file = "/Volumes/xchip_tcga/gbm/analysis/Barbara/Affy128_cnLOH/chr10.FALLS_p_TCGAaffxB4_1_GenomeWideSNP_6_A10_190536.cn.bin";
        dumpFloats(file);
    }

    public static void canonicalPath() {
        try {
            File f = new File("/Volumes/bob/sue/../../foo.txt");
            System.out.println(f.getCanonicalPath());
        } catch (IOException ex) {
            Logger.getLogger(MiscTests.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public static void printOS() {
        System.out.println(System.getProperty("os.name"));
        System.out.println(System.getProperty("os.arch"));

    }

    public static void compareFloats(File file1, File file2) throws FileNotFoundException, IOException {
        DataInputStream dis1 = new DataInputStream(new BufferedInputStream(new FileInputStream(file1)));
        DataInputStream dis2 = new DataInputStream(new BufferedInputStream(new FileInputStream(file2)));

        assert file1.length() == file2.length();

        int nFloats = (int) (file1.length() / 4);

        for (int i = 0; i < nFloats; i++) {
            float f1 = dis1.readFloat();
            float f2 = dis2.readFloat();

            if (f1 < 0) {
                if (!Float.isNaN(f2)) {
                    System.out.println("Difference! i=" + i + " f1= " + f1 + "  f2= " + f2);
                    return;
                }
            } else if (f1 != f2) {
                System.out.println("Difference! i=" + i + " f1= " + f1 + "  f2= " + f2);
                return;
            }
        }
    }

    public static void fixNaN(File inputDirectory, File outputDirectory) throws FileNotFoundException, IOException {
        if (!outputDirectory.exists()) {
            outputDirectory.mkdir();
        }

        for (File f : inputDirectory.listFiles()) {
            String fn = f.getName();

            if (fn.contains("snplocation") || fn.contains("snipid")) {

                // skip
            } else {
                DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(f)));
                int nFloats = (int) (f.length() / 4);
                File outputFile = new File(outputDirectory, f.getName());
                DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));

                for (int i = 0; i < nFloats; i++) {
                    float f1 = dis.readFloat();

                    if (f1 < 0) {
                        f1 = Float.NaN;
                    }

                    dos.writeFloat(f1);
                }

                dis.close();
                dos.close();
            }
        }
    }

    public static void dumpLongs(String file) throws FileNotFoundException, IOException {

        /*
         *  LocalBinaryReader reader = new LocalBinaryReader();
         * reader.setLittleEndian(true);
         * long[] longs = reader.readLongs(new File(file));
         * for (int i = 0; i < longs.length; i++) {
         * System.out.println(longs[i]);
         * }
         */
        DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(file)));

        for (int i = 0; i < 100; i++) {
            System.out.println(dis.readLong());
        }

        dis.close();
    }

    public static void dumpFloats(String file) throws FileNotFoundException, IOException {

        /*
         *  LocalBinaryReader reader = new LocalBinaryReader();
         * reader.setLittleEndian(true);
         * long[] longs = reader.readLongs(new File(file));
         * for (int i = 0; i < longs.length; i++) {
         * System.out.println(longs[i]);
         * }
         */
        DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(file)));

        for (int i = 0; i < 100; i++) {
            System.out.println(dis.readFloat());
        }

        dis.close();
    }


    /**
     * Print the last 4 floats
     *
     * @param file
     * @throws java.io.FileNotFoundException
     * @throws java.io.IOException
     */
    public static void tailLongs(String file) throws FileNotFoundException, IOException {

        /*
         *  LocalBinaryReader reader = new LocalBinaryReader();
         * reader.setLittleEndian(true);
         * long[] longs = reader.readLongs(new File(file));
         * for (int i = 0; i < longs.length; i++) {
         * System.out.println(longs[i]);
         * }
         */
        long len = (new File(file)).length();
        long pos = len - (8 * 32);

        DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(file)));
        dis.skip(pos);

        try {
            while (true) {
                System.out.println(dis.readLong());
            }

        } catch (EOFException eOFException) {
        }

        dis.close();
    }

    public static void moreLongs(String file) throws FileNotFoundException, IOException {

        /*
         *  LocalBinaryReader reader = new LocalBinaryReader();
         * reader.setLittleEndian(true);
         * long[] longs = reader.readLongs(new File(file));
         * for (int i = 0; i < longs.length; i++) {
         * System.out.println(longs[i]);
         * }
         */
        long len = (new File(file)).length();

        DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(file)));

        try {
            for (int i = 0; i < 10; i++) {
                System.out.println(dis.readLong());
            }

        } catch (EOFException eOFException) {
        }

        dis.close();
    }
}
