/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.sam.reader;

import htsjdk.samtools.*;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.sam.PicardAlignment;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
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
import java.util.zip.GZIPInputStream;

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

            SeekableStream indexStream = getInputStream(locator.getBamIndexPath());
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
                SeekableStream ss = new IGVSeekableBufferedStream(IGVSeekableStreamFactory.getInstance().getStreamFor(url));
                reader = getSamReader(locator, true);
            }
            CloseableIterator<SAMRecord> iter = reader.query(sequence, start + 1, end, contained);
            return new WrappedIterator(iter);
        } catch (IOException e) {
            log.error("Error opening SAM reader", e);
            throw new RuntimeException("Error opening SAM reader", e);
        }
    }


    private SeekableStream getInputStream(String indexPath) throws IOException {

        SeekableStream ss = null;
        URL indexURL = new URL(indexPath);
        boolean foundIndex;

        try {
            ss = IGVSeekableStreamFactory.getInstance().getStreamFor(indexURL);
            foundIndex = true;
        }
        catch(FileNotFoundException e) {

            String newIndexPath = indexPath.replace(".bam.bai", ".bai");
            indexURL = new URL(newIndexPath);
            try {
                ss = IGVSeekableStreamFactory.getInstance().getStreamFor(indexURL);
                foundIndex = true;
            }
            catch(FileNotFoundException e1) {
                foundIndex = false;
            }
        }


        if (!foundIndex) {
            String msg = "Index file not found: " + indexPath;
            throw new DataLoadException(msg, indexPath);
        }

        return ss;
    }


}
