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

package org.broad.igv.sam.reader;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.*;
import net.sf.samtools.util.SeekableBufferedStream;
import net.sf.samtools.util.SeekableStream;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.sam.Alignment;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.*;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.stream.IGVUrlHelper;
import org.broad.igv.util.stream.SeekablePicardStream;
import org.broad.tribble.util.SeekableFTPStream;

import java.io.*;
import java.net.URL;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Sep 22, 2009
 * Time: 2:21:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class BAMHttpQueryReader implements AlignmentQueryReader {

    static Logger log = Logger.getLogger(BAMHttpQueryReader.class);

    static Properties cachedEtags;

    URL url;
    SAMFileHeader header;
    File indexFile;
    SAMFileReader reader;

    public BAMHttpQueryReader(ResourceLocator locator, boolean requireIndex) throws IOException {
        this.url = new URL(locator.getPath());
        if (requireIndex) {
            indexFile = getIndexFile(url, locator.getIndexPath());
            if(indexFile == null) {
                throw new RuntimeException("Could not load index file for file: " + url.getPath());
            }
            SeekableStream ss = new SeekableBufferedStream(getSeekableStream(url));
            reader = new SAMFileReader(ss, indexFile, false);
        } else {
            InputStream is = HttpUtils.getInstance().openConnectionStream(url);
            reader = new SAMFileReader(new BufferedInputStream(is));
        }

    }

    public void close() throws IOException {
        if (reader != null) {
            reader.close();
        }
    }

    public SAMFileHeader getHeader() {
        if (header == null) {
            header = reader.getFileHeader();
        }
        return header;
    }

    public boolean hasIndex() {
        return indexFile != null && indexFile.exists();
    }

    public Set<String> getSequenceNames() {
        SAMFileHeader header = getHeader();
        if (header == null) {
            return null;
        }
        Set<String> seqNames = new HashSet();
        List<SAMSequenceRecord> records = header.getSequenceDictionary().getSequences();
        if (records.size() > 0) {
            for (SAMSequenceRecord rec : header.getSequenceDictionary().getSequences()) {
                String chr = rec.getSequenceName();
                seqNames.add(chr);
            }
        }
        return seqNames;
    }


    public CloseableIterator<Alignment> iterator() {
        try {
            if (reader == null) {
                InputStream is = HttpUtils.getInstance().openConnectionStream(url);
                reader = new SAMFileReader(new BufferedInputStream(is));
            }
            return new WrappedIterator(reader.iterator());
        } catch (IOException e) {
            log.error("Error creating iterator", e);
            throw new RuntimeException(e);
        }

    }

    public CloseableIterator<Alignment> query(String sequence, int start, int end, boolean contained) {
        try {
            if (reader == null) {
                SeekableStream ss = new SeekableBufferedStream(getSeekableStream(url));
                reader = new SAMFileReader(ss, indexFile, false);
            }
            CloseableIterator<SAMRecord> iter = reader.query(sequence, start, end, contained);
            return new WrappedIterator(iter);
        } catch (IOException e) {
            log.error("Error opening SAM reader", e);
            throw new RuntimeException("Error opening SAM reader", e);
        }
    }

    private SeekableStream getSeekableStream(URL url) throws IOException {
        String protocol = url.getProtocol().toLowerCase();
        SeekableStream is = null;
        if (protocol.equals("http") || protocol.equals("https")) {
            boolean useByteRange = HttpUtils.getInstance().useByteRange(url);
            if (useByteRange) {
                org.broad.tribble.util.SeekableStream tribbleStream =
                        new org.broad.tribble.util.SeekableHTTPStream(new IGVUrlHelper(url));
                String source = url.toExternalForm();
                is = new SeekablePicardStream(tribbleStream, source);
            } else {
                throw new RuntimeException("Byte-range requests are disabled.  HTTP and FTP access to BAM files require byte-range support.");
            }
        } else if (protocol.equals("ftp")) {
            org.broad.tribble.util.SeekableStream tribbleStream = new SeekableFTPStream(url);
            String source = url.toExternalForm();
            is = new SeekablePicardStream(tribbleStream, source);
        } else {
            throw new RuntimeException("Unknown protocol: " + protocol);
        }
        return is;
    }


    // TODO -- revisit caching scehme,  do something for ftp loads

    File getIndexFile(URL url, String indexPath) throws IOException {

        String urlString = url.toString();

        // Create a filename unique for this url;
        String idxFilename = getTmpIndexFilename(urlString);

        indexFile = new File(this.getCacheDirectory(), idxFilename);
        if (indexFile.exists()) {
            indexFile.delete();
        }
        loadIndexFile(urlString, indexPath, indexFile);
        indexFile.deleteOnExit();

        return indexFile;

    }

    private String getTmpIndexFilename(String bamURL) {
        int tmp = bamURL.lastIndexOf("/");
        String prefix = tmp > 0 ? bamURL.substring(tmp + 1) : "index_";
        String indexName = prefix + System.currentTimeMillis() + ".bai";
        return indexName;
    }

    private void loadIndexFile(String path, String indexPath, File indexFile) throws IOException {
        InputStream is = null;
        OutputStream os = null;

        try {
            String idx = (indexPath != null && indexPath.length() > 0) ? indexPath : path + ".bai";
            URL indexURL = new URL(idx);
            os = new FileOutputStream(indexFile);
            try {
                is = HttpUtils.getInstance().openConnectionStream(indexURL);
            } catch (FileNotFoundException e) {
                // Try other index convention
                String baseName = path.substring(0, path.length() - 4);
                indexURL = new URL(baseName + ".bai");

                try {
                    is = org.broad.igv.util.HttpUtils.getInstance().openConnectionStream(indexURL);
                } catch (FileNotFoundException e1) {
                    MessageUtils.showMessage("Index file not found for file: " + path);
                    throw new DataLoadException("Index file not found for file: " + path, path);
                }
            }
            byte[] buf = new byte[512000];
            int bytesRead;
            while ((bytesRead = is.read(buf)) != -1) {
                os.write(buf, 0, bytesRead);
            }

        } finally {
            if (is != null) {
                try {
                    is.close();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }
            if (os != null) {
                try {
                    os.close();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }

        }
    }


    // TODO -- move everything below to a utility class
    static File cacheDirectory = null;

    static synchronized File getCacheDirectory() {
        if (cacheDirectory == null) {
            File defaultDir = Globals.getIgvDirectory();
            if (defaultDir.exists()) {
                cacheDirectory = new File(defaultDir, "bam");
                if (!cacheDirectory.exists()) {
                    cacheDirectory.mkdir();
                }
            }
        }
        return cacheDirectory;
    }

    static Properties tagCache = null;
    static Properties indexLookup = null;

    static synchronized Properties getTagCache() {
        if (tagCache == null) {
            tagCache = loadCachedProperties("etags");
        }
        return tagCache;
    }

    static synchronized Properties getIndexCache() {
        if (indexLookup == null) {
            indexLookup = loadCachedProperties("index");
        }
        return indexLookup;
    }

    private static Properties loadCachedProperties(String name) {
        File propFile = new File(getCacheDirectory(), name);
        Properties props = new Properties();
        if (propFile.exists()) {
            InputStream reader = null;
            try {
                reader = new FileInputStream(propFile);
                props.load(reader);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            } finally {
                if (reader != null) {
                    try {
                        reader.close();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
        }
        return props;
    }


    private static void writeCachedProperties(Properties props, String name) {
        OutputStream os = null;
        try {
            File propFile = new File(getCacheDirectory(), name);
            os = new FileOutputStream(propFile);
            props.store(os, null);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } finally {
            if (os != null) {
                try {
                    os.close();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }
        }
    }
}
