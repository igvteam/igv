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
