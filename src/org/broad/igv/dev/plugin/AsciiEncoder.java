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

import org.broad.tribble.Feature;

import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * User: jacob
 * Date: 2012-Oct-01
 */
public class AsciiEncoder<E extends Feature> implements FeatureEncoder<E> {

    protected LineFeatureEncoder<E> lineFeatureEncoder;

    public AsciiEncoder(LineFeatureEncoder<E> lineFeatureEncoder) {
        this.lineFeatureEncoder = lineFeatureEncoder;
    }

    public static final String NUM_COLS_ATTR = "numCols";

    @Override
    public Map<String, Object> encodeAll(OutputStream outputStream, Iterator<E> features) {
        LineFeatureEncoder<E> encoder = lineFeatureEncoder;
        PrintWriter writer = new PrintWriter(new OutputStreamWriter(outputStream));
        Map<String, Object> attributes = new HashMap<String, Object>(1);
        int tmpNumCols = 0;

        if (features != null) {
            String header = encoder.getHeader();
            if (header != null) {
                writer.println(header);
            }
            E feature;
            while (features.hasNext()) {
                feature = features.next();
                String line = encoder.encode(feature);
                if (line == null) continue;
                writer.println(line);

                //We do not require consistency of output
                tmpNumCols = encoder.getNumCols(line);
            }
        }

        attributes.put(NUM_COLS_ATTR, tmpNumCols);
        writer.flush();
        writer.close();
        return attributes;
    }

    @Override
    public void setInputs(List<String> commands, Map<Argument, Object> argumentMap) {
        lineFeatureEncoder.setInputs(commands, argumentMap);
    }
}
