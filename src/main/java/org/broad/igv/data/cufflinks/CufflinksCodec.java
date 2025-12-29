package org.broad.igv.data.cufflinks;

import org.broad.igv.logging.*;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.util.ParsingUtils;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.readers.LineIterator;

/**
 * @author jacob
 * @date 2013-Apr-18
 */
public abstract class CufflinksCodec<T extends Feature> extends AsciiFeatureCodec<T> {

    private static Logger log = LogManager.getLogger(CufflinksCodec.class);

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

