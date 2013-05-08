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
import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


/**
 * @author jrobinso
 * @date Jun 24, 2010
 */
public class CachingFeatureReader extends AbstractCacher implements IGVFeatureReader {

    private static Logger log = Logger.getLogger(CachingFeatureReader.class);
    private static int maxBinCount = 1000;
    private static int defaultBinSize = 16000; // <= 16 kb

    private FeatureReader tribbleFeatureReader;


    public CachingFeatureReader(FeatureReader tribbleFeatureReader) {
        this(tribbleFeatureReader, maxBinCount, defaultBinSize);
    }


    public CachingFeatureReader(FeatureReader tribbleFeatureReader, int binCount, int binSize) {
        super(binCount, binSize);
        this.tribbleFeatureReader = tribbleFeatureReader;
    }


    @Override
    protected Iterator<Feature> queryRaw(String chr, int start, int end) throws IOException {

        // Tribble iterators must be closed, so we need to copy the features and insure closure before exiting.
        CloseableTribbleIterator<Feature> iter = null;
        try {
            iter = tribbleFeatureReader.query(chr, start, end);
            List<Feature> featureList = new ArrayList<Feature>();
            while (iter.hasNext()) {
                Feature f = iter.next();
                if (f.getStart() > end) {
                    break;
                } else if (f.getEnd() < start) {
                    continue;
                } else {
                    featureList.add(f);
                }
            }
            return featureList.iterator();
        } finally {
            if (iter != null) iter.close();
        }
    }


    @Override
    public List<String> getSequenceNames() {
        return tribbleFeatureReader.getSequenceNames();
    }

    @Override
    public Object getHeader() {
        return tribbleFeatureReader.getHeader();
    }

    @Override
    public Iterator query(String chr, int start, int end) throws IOException {
        return queryCached(chr, start, end);
    }

    /**
     * Return an iterator over the entire file.  Nothing to cache,  just delegate to the wrapped reader
     *
     * @throws java.io.IOException
     */
    @Override
    public Iterator iterator() throws IOException {
        return tribbleFeatureReader.iterator();
    }


}

