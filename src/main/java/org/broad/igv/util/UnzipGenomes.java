/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.util;


import java.io.*;
import java.util.*;
import java.util.zip.*;
import java.text.*;


/**
 * @author jrobinso
 * @date Dec 22, 2010
 */
public class UnzipGenomes {


    public static void main(String[] args) throws IOException {
        File root = new File("/Users/jrobinso/projects/genomes");

        for (File f : root.listFiles()) {

            if (f.getName().endsWith(".genome")) {
                String dirName = f.getName().replace(".genome", "");
                File dir = new File(root, dirName);
                dir.mkdir();
                unzip(f, dir);
            }
        }
    }


    public static void unzip(File zipFile, File destDir) throws IOException {
        InputStream in = new BufferedInputStream(new FileInputStream(zipFile));
        ZipInputStream zin = new ZipInputStream(in);
        ZipEntry e;
        while ((e = zin.getNextEntry()) != null) {

            File f = new File(destDir, e.getName());
            unzip(zin, f);
        }
        zin.close();
    }

    public static void unzip(ZipInputStream zin, File f)
            throws IOException {
        System.out.println("unzipping " + f.getAbsolutePath());
        FileOutputStream out = null;
        try {
            out = new FileOutputStream(f);
            byte[] b = new byte[512];
            int len = 0;
            while ((len = zin.read(b)) != -1) {
                out.write(b, 0, len);
            }
        } catch (FileNotFoundException e) {
            System.err.println("Error unzipping: " + f.getAbsolutePath());
        }
        finally {
            if (out != null)
                out.close();
        }
    }

}
