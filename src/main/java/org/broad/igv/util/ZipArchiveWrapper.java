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
