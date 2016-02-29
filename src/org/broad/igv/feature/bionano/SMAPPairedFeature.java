/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 University of California San Diego
 * Author: Jim Robinson
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

package org.broad.igv.feature.bionano;

import org.apache.log4j.Logger;
import org.broad.igv.feature.AbstractFeature;
import org.broad.igv.track.WindowFunction;


/**
 * Created by jrobinson on 2/25/16.
 */
public class SMAPPairedFeature extends AbstractFeature {

    private static Logger log = Logger.getLogger(SMAPPairedFeature.class);

    SMAPFeature feature1;
    SMAPFeature feature2;

    String type;

    public SMAPPairedFeature(SMAPFeature feature1, SMAPFeature feature2) {


        if (!feature1.getChr().equals(feature2.getChr())) {
            // TODO - throw error?
            log.error("Inter-chromosomal linked features not supported");
            return;
        }

        setChr(feature1.getChr());
        if (feature1.getStart() < feature2.getStart()) {
            this.feature1 = feature1;
            this.feature2  = feature2;
        } else {
            this.feature1 = feature2;
            this.feature2  = feature1;
        }

        // Feaures shouldn't overlap, but just in case check both
        setStart(Math.min(feature1.getStart(), feature2.getStart()));
        setEnd(Math.max(feature1.getEnd(), feature2.getEnd()));

    }


    @Override
    public String getValueString(double position, WindowFunction windowFunction) {

        StringBuffer buf = new StringBuffer();
        buf.append(feature1.getValueString(position, windowFunction));
        buf.append("<hr>");
        buf.append(feature2.getValueString(position, windowFunction));
        return buf.toString();

    }

    @Override
    public String getURL() {
        return null;
    }
}

