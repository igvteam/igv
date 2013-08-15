/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.data.cufflinks;

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.readers.LineIterator;

/**
 * @author jacob
 * @date 2013-Apr-18
 */
public abstract class CufflinksCodec<T extends Feature> extends AsciiFeatureCodec<T> {

    private static Logger log = Logger.getLogger(CufflinksCodec.class);

    String path;

    protected CufflinksCodec(Class<T> clazz, String path){
        super(clazz);
        this.path = path;
    }

    protected abstract Object readHeader(String[] tokens);

    @Override
    public Object readActualHeader(LineIterator reader){
        String headerLine = null;
        try {
            headerLine = reader.next();
            String[] tokens = ParsingUtils.TAB_PATTERN.split(headerLine);
            return readHeader(tokens);
        } catch (Exception e) {
            log.error(e.getMessage(), e);
            throw new DataLoadException("Error reading header: " + e.getMessage(), this.path);
        }
    }
}

