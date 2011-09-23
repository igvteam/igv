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

import org.apache.log4j.Logger;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;
import java.util.Iterator;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;


/**
 * @author eflakes
 */

final public class ZipArchiveWrapper {

    private static Logger logger = Logger.getLogger(Utilities.class);
    private File zipFile;
    private URL zipUrl;

    public ZipArchiveWrapper(File zipFile) throws FileNotFoundException, IOException {

        if (zipFile == null) {
            throw new FileNotFoundException("Zip file was null!");
        }
        this.zipFile = zipFile;
        zipUrl = zipFile.toURI().toURL();
    }

    public ZipIterator iterator() throws IOException {
        ZipIterator iterator =
                new ZipIterator(new ZipInputStream(HttpUtils.getInstance().openConnectionStream(zipUrl)));
        return iterator;
    }

    public int getEntryCount() throws IOException {
        ZipIterator iterator =
                new ZipIterator(new ZipInputStream(HttpUtils.getInstance().openConnectionStream(zipUrl)));
        int count = 0;
        while (iterator.hasNext()) {
            count++;
        }
        return count;
    }
    /*
     * Iterator
     */
    /*
    public class ZipIterator implements Iterator {

        public boolean hasNext() {
            if(zipInputStream == null) {
                return false;
            }
            else {
                try {
                    return zipInputStream.available() == 1 ? true:false;
                } catch(IOException e) {
                    logger.error("Zip entry error!", e);
                    return false;
                }
            }
        }

        public ZipEntry next() {

            ZipEntry entry = null;
            try {     
                entry = zipInputStream.getNextEntry();                         
            }
            catch(IOException e) {
                logger.error("Zip entry error!", e);
            }
            return entry;           
        }
*/


    public class ZipIterator implements Iterator {

        private boolean firstTime = true;
        private ZipEntry zipEntry = null;
        ZipInputStream zipInputStream;

        public ZipIterator(ZipInputStream zipInputStream) {
            this.zipInputStream = zipInputStream;
        }

        public boolean hasNext() {

            if (zipInputStream == null) {
                return false;
            } else {
                try {
                    if (firstTime) {
                        firstTime = false;
                        zipEntry = zipInputStream.getNextEntry();
                        if (zipEntry == null) {
                            return false; // must have read the last one already
                        }
                    } else {
                        if (zipEntry == null) {
                            return false; // must have read the last one already
                        } else {
                            zipEntry = zipInputStream.getNextEntry();
                            if (zipEntry == null) {
                                return false; // No more entries
                            }
                        }
                    }
                } catch (IOException e) {
                    logger.error("Zip entry error!", e);
                }
            }
            return true;
        }

        public ZipEntry next() {
            return zipEntry;
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }

        public ZipInputStream getZipInputStream() {
            return zipInputStream;
        }

        public void close() {
            try {
                if (zipInputStream != null) {
                    zipInputStream.close();
                }
            } catch (IOException ex) {
                logger.warn("Error closing zip file " + zipFile.getAbsolutePath(), ex);
            }
        }
    }

}
