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
import java.util.Iterator;

/**
 * User: jacob
 * Date: 2012-Oct-01
 */
public class AsciiEncoder<E extends Feature> implements FeatureEncoder<E> {

    protected LineFeatureEncoder<E> lineFeatureEncoder;

    public AsciiEncoder(LineFeatureEncoder<E> lineFeatureEncoder) {
        this.lineFeatureEncoder = lineFeatureEncoder;
    }

    @Override
    public int encodeAll(OutputStream outputStream, Iterator<E> features) {
        LineFeatureEncoder<E> encoder = lineFeatureEncoder;
        PrintWriter writer = new PrintWriter(new OutputStreamWriter(outputStream));

        int allNumCols = -1;
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

                //We require consistency of output
                int tmpNumCols = encoder.getNumCols(line);
                if (allNumCols < 0) {
                    allNumCols = tmpNumCols;
                } else {
                    assert tmpNumCols == allNumCols;
                }
            }
        }
        writer.flush();
        writer.close();
        return allNumCols;
    }
}
