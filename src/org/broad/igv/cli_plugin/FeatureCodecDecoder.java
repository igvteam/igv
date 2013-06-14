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

package org.broad.igv.cli_plugin;

import org.apache.log4j.Logger;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.readers.PositionalBufferedStream;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * @author jacob
 * @date 2013-Jun-14
 */
public class FeatureCodecDecoder<D extends Feature> implements FeatureDecoder<D> {

    private static Logger log = Logger.getLogger(PluginSource.class);

    private FeatureCodec<D> codec;

    public FeatureCodecDecoder(FeatureCodec<D> codec) {
        this.codec = codec;
    }

    @Override
    public Iterator<D> decodeAll(InputStream is, boolean strictParsing) throws IOException {
        PositionalBufferedStream pis = new PositionalBufferedStream(is);

        this.codec.readHeader(pis);

        List<D> features = new ArrayList<D>();
        while (!pis.isDone()) {
            try {
                D feature = this.codec.decode(pis);
                if (feature != null) {
                    features.add(feature);
                }
            } catch (Exception e) {
                log.error(e.getMessage(), e);
                if (strictParsing) {
                    throw new RuntimeException(e);
                }
            } finally {
                pis.close();
            }
        }
        return features.iterator();
    }

    @Override
    public void setAttributes(List<Map<String, Object>> attributes) {
        //pass
    }

    @Override
    public void setInputs(List<String> commands, Map<Argument, Object> argumentMap) {
        //pass
    }

}
