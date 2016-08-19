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

package org.broad.igv.feature;

import org.broad.igv.track.WindowFunction;

import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 7/18/13
 *         Time: 9:50 PM
 */
public class GFFFeature extends BasicFeature {

    List<String> componentAttributes = new ArrayList<>();

    public GFFFeature(BasicFeature feature) {
        super(feature);
    }

    public GFFFeature(String chr, int start, int end, Strand strand) {
        super(chr, start, end, strand);
    }

    @Override
    public void addExon(Exon region) {
        super.addExon(region);
    }

    @Override
    public void addUTRorCDS(BasicFeature bf) {
        super.addUTRorCDS(bf);
    }

    @Override
    public String getValueString(double position, int mouseX, WindowFunction ignored) {

        StringBuffer valueString = new StringBuffer();

        valueString.append("<b>Type:</b>&nbsp;" + this.type);

        if (attributes != null) {
            valueString.append(getAttributeString());
        }

        if (componentAttributes.size() > 0) {
            for (String s : componentAttributes) {
                valueString.append("---------------------------");
                valueString.append(s);
            }
        }
        return valueString.toString();
    }

    public void mergeAttributes(BasicFeature mrnaPart) {
        StringBuffer buf = new StringBuffer();
        buf.append("<br><b>Type:</b>&nbsp;" + mrnaPart.getType());
        mrnaPart.getAttributes().printHtml(buf, 100);
        componentAttributes.add(buf.toString());
    }
}
