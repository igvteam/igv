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

package org.broad.igv.sam.reader;

import htsjdk.samtools.*;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.sam.PicardAlignment;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.stream.IGVSeekableBufferedStream;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Sep 22, 2009
 * Time: 2:21:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class BAMHttpReader implements AlignmentReader<PicardAlignment> {

    static Logger log = Logger.getLogger(BAMHttpReader.class);

    // Length of day in milliseconds
    public static final long oneDay = 24 * 60 * 60 * 1000;

    static Hashtable<String, File> indexFileCache = new Hashtable<String, File>();
    private final ResourceLocator locator;

    URL url;
    SAMFileHeader header;
    htsjdk.samtools.SamReader reader;
    List<String> sequenceNames;
    private boolean indexed = false; // False until proven otherwise

    public BAMHttpReader(ResourceLocator locator, boolean requireIndex) throws IOException {
        this.locator = locator;
        this.url = new URL(locator.getPath());
        reader = getSamReader(locator, requireIndex);

    }

    public SamReader getSamReader(ResourceLocator locator, boolean requireIndex) throws IOException {

        if (requireIndex) {
            final SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);

            SeekableStream indexStream = getIndexStream(locator.getBamIndexPath());
            this.indexed = true;

            SeekableStream ss = new IGVSeekableBufferedStream(IGVSeekableStreamFactory.getInstance().getStreamFor(url), 128000);
            SamInputResource resource = SamInputResource.of(ss).index(indexStream);
            return factory.open(resource);
        } else {
            InputStream is = HttpUtils.getInstance().openConnectionStream(url);
            return new SAMFileReader(new BufferedInputStream(is));
        }
    }

    public void close() throws IOException {
        if (reader != null) {
            reader.close();
        }
    }

    public SAMFileHeader getFileHeader() {
        if (header == null) {
            header = reader.getFileHeader();
        }
        return header;
    }

    public boolean hasIndex() {
        return indexed;
    }

    public Set<String> getPlatforms() {
        return AlignmentReaderFactory.getPlatforms(getFileHeader());
    }


    public List<String> getSequenceNames() {
        if (sequenceNames == null) {
            SAMFileHeader header = getFileHeader();
            if (header == null) {
                return null;
            }
            sequenceNames = new ArrayList();
            List<SAMSequenceRecord> records = header.getSequenceDictionary().getSequences();
            if (records.size() > 0) {
                for (SAMSequenceRecord rec : header.getSequenceDictionary().getSequences()) {
                    String chr = rec.getSequenceName();
                    sequenceNames.add(chr);
                }
            }
        }
        return sequenceNames;
    }


    public CloseableIterator<PicardAlignment> iterator() {
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

    public CloseableIterator<PicardAlignment> query(String sequence, int start, int end, boolean contained) {
        try {
            if (reader == null) {
                reader = getSamReader(locator, true);
            }
            CloseableIterator<SAMRecord> iter = reader.query(sequence, start + 1, end, contained);
            return new WrappedIterator(iter);
        } catch (IOException e) {
            log.error("Error opening SAM reader", e);
            throw new RuntimeException("Error opening SAM reader", e);
        }
    }


    private SeekableStream getIndexStream(String indexPath) throws IOException {

        List<String> pathsTried = new ArrayList<String>();

        // Strip parameters


        if (HttpUtils.getInstance().resourceAvailable(new URL(indexPath))) {
            return IGVSeekableStreamFactory.getInstance().getStreamFor(new URL(indexPath));
        } else {


            if (indexPath.contains(".bam.bai")) {
                indexPath = indexPath.replaceFirst(".bam.bai", ".bai");
                pathsTried.add(indexPath);
                if (HttpUtils.getInstance().resourceAvailable(new URL(indexPath))) {
                    log.info("Index found: " + indexPath);
                    return IGVSeekableStreamFactory.getInstance().getStreamFor(new URL(indexPath));
                }
            } else {
                indexPath = indexPath.replaceFirst(".bai", ".bam.bai");
                pathsTried.add(indexPath);
                if (HttpUtils.getInstance().resourceAvailable(new URL(indexPath))) {
                    log.info("Index found: " + indexPath);
                    return IGVSeekableStreamFactory.getInstance().getStreamFor(new URL(indexPath));
                }
            }
        }

        String msg = "Index file not found.  Tried ";
        for (String path : pathsTried) {
            msg += "<br>" + indexPath;
        }
        throw new DataLoadException(msg, indexPath);

    }


}
