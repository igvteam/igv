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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.data;

import org.broad.igv.util.FileUtils;

import java.io.*;
import java.util.zip.CRC32;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

/**
 * @author jrobinso
 */
public class ZipUtils {

    public static void zipDirectory(File dir2Zip, File outputFile) {
        try {
            //create a ZipOutputStream to zip the data to 
            ZipOutputStream zos = new ZipOutputStream(new FileOutputStream(outputFile));
            //assuming that there is a directory named inFolder (If there 
            //isn't create one) in the same directory as the one the code 

            //call the zipDir method 
            zipDirRecursive(dir2Zip, dir2Zip, zos);
            //close the stream 
            zos.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    //here is the code for the method
    public static void zipDirRecursive(File zipDir, File currentDir, ZipOutputStream zos) {
        try {
            //create a new File object based on the directory we have to zip 

            //get a listing of the directory content 
            String[] dirList = zipDir.list();
            byte[] readBuffer = new byte[2156];
            int bytesIn = 0;
            //loop through dirList, and zip the files 
            for (File f : currentDir.listFiles()) {
                if (f.isDirectory()) {
                    //if the File object is a directory, call this 
                    //function again to add its content recursively 
                    zipDirRecursive(zipDir, f, zos);
                    //loop again 
                    continue;
                } else {
                    //if we reached here, the File object f was not a directory 
                    //create a FileInputStream on top of f 

                    FileInputStream fis = new FileInputStream(f);

                    //create a  new zipentry.  Use a relative path
                    String relativePath = FileUtils.getRelativePath(zipDir.getAbsolutePath(), f.getAbsolutePath(),
                            System.getProperty("file.separator"));
                    String piPath = FileUtils.getPlatformIndependentPath(relativePath);
                    if (piPath.startsWith("./")) {
                        piPath = piPath.substring(2);
                    }

                    ZipEntry anEntry = new ZipEntry(piPath);
                    anEntry.setMethod(ZipEntry.STORED);
                    anEntry.setCompressedSize(f.length());
                    anEntry.setSize(f.length());
                    anEntry.setCrc(getCrc(f));

                    zos.putNextEntry(anEntry);
                    //now write the content of the file to the ZipOutputStream 
                    while ((bytesIn = fis.read(readBuffer)) != -1) {
                        zos.write(readBuffer, 0, bytesIn);
                    }
                    //close the Stream 
                    fis.close();
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static long getCrc(File file) throws IOException {
        CRC32 crc = new CRC32();
        crc.reset();
        byte[] buffer = new byte[1024];
        int bytesRead = 0;
        BufferedInputStream bis = new BufferedInputStream(new FileInputStream(file));

        while ((bytesRead = bis.read(buffer)) != -1) {
            crc.update(buffer, 0, bytesRead);
        }
        bis.close();
        return crc.getValue();

    }
}
