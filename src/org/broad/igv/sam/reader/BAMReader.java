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
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.stream.IGVSeekableBufferedStream;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Sep 22, 2009
 * Time: 2:21:04 PM
 */
public class BAMReader implements AlignmentReader<PicardAlignment> {

    static Logger log = Logger.getLogger(BAMReader.class);

    private final ResourceLocator locator;

    SAMFileHeader header;
    htsjdk.samtools.SamReader reader;
    List<String> sequenceNames;
    private boolean indexed = false; // False until proven otherwise

    public BAMReader(ResourceLocator locator, boolean requireIndex) throws IOException {
        this.locator = locator;
        reader = getSamReader(locator, requireIndex);
        header = reader.getFileHeader();
    }

    private SamReader getSamReader(ResourceLocator locator, boolean requireIndex) throws IOException {

        boolean isLocal = locator.isLocal();
        final SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        SamInputResource resource;

        if (isLocal) {
            resource = SamInputResource.of(new File(locator.getPath()));
        } else {
            URL url = new URL(locator.getPath());
            if(requireIndex) {
                resource = SamInputResource.of(new IGVSeekableBufferedStream(IGVSeekableStreamFactory.getInstance().getStreamFor(url), 128000));
            }
            else {
                resource = SamInputResource.of(HttpUtils.getInstance().openConnectionStream(url));
            }
        }

        if (requireIndex) {

            String indexPath = getIndexPath(locator);
            indexed = true;
            if (isLocal) {
                File indexFile = new File(indexPath);
                resource = resource.index(indexFile);
            } else {
                SeekableStream indexStream = IGVSeekableStreamFactory.getInstance().getStreamFor(new URL(indexPath));
                resource = resource.index(indexStream);
            }
        }

        return factory.open(resource);

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
        return new WrappedIterator(reader.iterator());
    }

    public CloseableIterator<PicardAlignment> query(String sequence, int start, int end, boolean contained) {
        CloseableIterator<SAMRecord> iter = reader.query(sequence, start + 1, end, contained);
        return new WrappedIterator(iter);
    }

    private String getIndexPath(ResourceLocator locator) throws IOException {

        List<String> pathsTried = new ArrayList<String>();

        String path = locator.getPath();
        String indexPath = locator.getIndexPath();


        if (indexPath != null) {
            return indexPath;  // Explicit
        } else {

            if (path.toLowerCase().startsWith("http://") || path.toLowerCase().startsWith("https://")) {
                // See if bam file is specified by parameter
                try {
                    URL url = new URL(path);
                    String queryString = url.getQuery();
                    if (queryString != null) {
                        Map<String, String> parameters = HttpUtils.parseQueryString(queryString);
                        if (parameters.containsKey("index")) {
                            indexPath = parameters.get("index");
                        } else if (parameters.containsKey("file")) {
                            String bamFile = parameters.get("file");
                            String bamIndexFile = bamFile + ".bai";
                            String newQueryString = queryString.replace(bamFile, bamIndexFile);
                            indexPath = path.replace(queryString, newQueryString);
                        } else {
                            indexPath = path.replace(url.getPath(), url.getPath() + ".bai");
                        }
                    }
                } catch (MalformedURLException e) {
                    log.error(e.getMessage(), e);
                }
            }
        }
        if (indexPath != null && FileUtils.resourceExists(indexPath)) {
            return indexPath;
        }

        if (indexPath == null) {
            indexPath = path + ".bai";
        }

        if (FileUtils.resourceExists(indexPath)) {
            return indexPath;
        } else {
            if (indexPath.contains(".bam.bai")) {
                indexPath = indexPath.replaceFirst(".bam.bai", ".bai");
                pathsTried.add(indexPath);
                if (FileUtils.resourceExists(indexPath)) {
                    return indexPath;
                }
            } else {
                indexPath = indexPath.replaceFirst(".bai", ".bam.bai");
                pathsTried.add(indexPath);
                if (FileUtils.resourceExists(indexPath)) {
                    return indexPath;
                }
            }
        }

        if (indexPath == null) {
            String defaultValue = locator.getPath() + ".bai";
            indexPath = MessageUtils.showInputDialog(
                    "Index is required, but no index found.  Please enter path to index file:",
                    defaultValue);
            if (indexPath != null && FileUtils.resourceExists(indexPath)) {
                return indexPath;
            }
        }

        String msg = "Index file not found.  Tried ";
        for (String p : pathsTried) {
            msg += "<br>" + p;
        }
        throw new DataLoadException(msg, indexPath);

    }


}
