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
