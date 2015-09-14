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

package org.broad.igv.cli_plugin;

import htsjdk.tribble.Feature;

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

    /**
     * Current use-case is BEDTools, which requires unix-style
     * line endings (linefeed). Instead of ending the line
     * based on platform, we use a constant.
     */
    public static final String EOL_CHAR = "\n";

    protected LineFeatureEncoder<E> lineFeatureEncoder;

    public AsciiEncoder(LineFeatureEncoder<E> lineFeatureEncoder) {
        this.lineFeatureEncoder = lineFeatureEncoder;
    }

    public static final String NUM_COLS_ATTR = "numCols";

    @Override
    public Map<String, Object> encodeAll(OutputStream outputStream, Iterator<? extends E> features) {
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
                writer.print(line);
                writer.print(EOL_CHAR);

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
    public void setInputs(List<String> commands, Map<Argument, Object> argumentMap, Argument argument) {
        lineFeatureEncoder.setInputs(commands, argumentMap, argument);
    }
}
