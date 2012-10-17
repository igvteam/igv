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

package org.broad.igv.dev.plugin;

import org.apache.log4j.Logger;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * User: jacob
 * Date: 2012-Sep-27
 */
public class AsciiDecoder<D extends Feature> implements FeatureDecoder<D> {

    private static Logger log = Logger.getLogger(AsciiDecoder.class);

    protected LineFeatureDecoder<D> lineFeatureDecoder;

    public AsciiDecoder() {
    }

    public AsciiDecoder(LineFeatureDecoder<D> lineFeatureDecoder) {
        this.lineFeatureDecoder = lineFeatureDecoder;
    }

    public Iterator<D> decodeAll(InputStream is, boolean strictParsing) throws IOException {

        List<D> featuresList = new ArrayList<D>();
        BufferedReader bis = new BufferedReader(new InputStreamReader(is));
        String line;
        D feat;

        while ((line = bis.readLine()) != null) {
            try {
                feat = decode(line);
                if (feat != null) {
                    featuresList.add(feat);
                }
            } catch (Exception e) {
                log.error(e);
                if (strictParsing) {
                    throw new RuntimeException(e);
                }
            }
        }

        is.close();
        return featuresList.iterator();
    }

    public D decode(String line) {
        return this.lineFeatureDecoder.decode(line);
    }

    @Override
    public void setAttributes(List<Map<String, Object>> attributes) {
    }

    @Override
    public void setInputs(List<String> commands, Map<Argument, Object> argumentMap) {
    }

    /**
     * Wrap an AsciiFeatureCodec into implementing LineFeatureDecoder
     *
     * @param <T>
     */
    public static class DecoderWrapper<T extends Feature> extends AsciiDecoder<T> implements LineFeatureDecoder<T> {

        private AsciiFeatureCodec<T> wrappedCodec;

        public DecoderWrapper(AsciiFeatureCodec<T> wrappedCodec) {
            this.wrappedCodec = wrappedCodec;
            this.lineFeatureDecoder = this;
        }

        @Override
        public T decode(String line) {
            return wrappedCodec.decode(line);
        }

    }
}
