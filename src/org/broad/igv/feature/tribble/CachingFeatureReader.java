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


package org.broad.igv.feature.tribble;

import org.apache.log4j.Logger;
import org.broad.igv.feature.AbstractCacher;
import org.broad.igv.feature.WrappedIterator;
import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureReader;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;


/**
 * @author jrobinso
 * @date Jun 24, 2010
 */
public class CachingFeatureReader extends AbstractCacher implements FeatureReader {

    private static Logger log = Logger.getLogger(CachingFeatureReader.class);
    private static int maxBinCount = 1000;
    private static int defaultBinSize = 16000; // <= 16 kb

    private FeatureReader reader;


    public CachingFeatureReader(FeatureReader reader) {
        this(reader, maxBinCount, defaultBinSize);
    }


    public CachingFeatureReader(FeatureReader reader, int binCount, int binSize) {
        super(binCount, binSize);
        this.reader = reader;
    }


    @Override
    protected Iterator<Feature> queryRaw(String chr, int start, int end) throws IOException {
        return reader.query(chr, start, end);
    }

    @Override
    public List<String> getSequenceNames() {
        return reader.getSequenceNames();
    }

    @Override
    public Object getHeader() {
        return reader.getHeader();
    }

    @Override
    public CloseableTribbleIterator query(String chr, int start, int end) throws IOException {
        return new WrappedIterator(queryCached(chr, start, end));
    }

    /**
     * Return an iterator over the entire file.  Nothing to cache,  just delegate to the wrapped reader
     *
     * @throws java.io.IOException
     */
    @Override
    public CloseableTribbleIterator iterator() throws IOException {
        return reader.iterator();
    }

    @Override
    public void close() throws IOException {
        super.close();
        reader.close();
    }

}

