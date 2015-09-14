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

package org.broad.igv.feature.tribble;

import org.broad.igv.renderer.SpliceJunctionRenderer;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ParsingUtils;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.exception.CodecLineParsingException;
import htsjdk.tribble.readers.LineIterator;

/**
 * @author jrobinso
 * @date Aug 5, 2010
 */
public abstract class UCSCCodec<T extends Feature> extends AsciiFeatureCodec<T> {

    GFFCodec.GFF3Helper tagHelper = new GFFCodec.GFF3Helper();
    protected boolean gffTags = false;
    protected boolean spliceJunctions;

    FeatureFileHeader header;

    protected UCSCCodec(Class myClass) {
        super(myClass);
    }

    /**
     * @param reader
     * @return
     */
    @Override
    public Object readActualHeader(LineIterator reader) {
        String line;
        try {
            while (reader.hasNext()) {
                line = reader.peek();
                if (line.startsWith("#") || line.startsWith("track") ||
                        line.startsWith("browser")) {
                    readHeaderLine(line);
                    reader.next();
                } else {
                    break;
                }
            }
            return header;
        } catch (Exception e) {
            throw new CodecLineParsingException("Error parsing header", e);
        }
    }

    /**
     * Extract information from the header line.
     * Side effects: Calling this will create a new header field
     * if one is null. In general, should check whether the line
     * is a header line or not first.
     *
     * @param line
     * @return True iff any information was retrieved.
     */
    protected boolean readHeaderLine(String line) {
        //Header line found, may not have any content
        if (header == null) {
            header = new FeatureFileHeader();
        }
        if (line.startsWith("#type")) {
            String[] tokens = line.split("=");
            if (tokens.length > 1) {
                try {
                    header.setTrackType(TrackType.valueOf(tokens[1]));
                } catch (Exception e) {
                    // log.error("Error converting track type: " + tokens[1]);
                }
            }
        } else if (line.startsWith("#track") || line.startsWith("track")) {
            TrackProperties tp = new TrackProperties();
            ParsingUtils.parseTrackLine(line, tp);
            header.setTrackProperties(tp);
            if(gffTags == false) gffTags = tp.isGffTags();

            Class rendererClass = tp.getRendererClass();
            if (rendererClass != null && rendererClass.isAssignableFrom(SpliceJunctionRenderer.class)) {
                spliceJunctions = true;
            }

        } else if (line.toLowerCase().contains("#gfftags")) {
            gffTags = true;
        } else {
            return false;
        }
        return true;
    }

    public Feature decodeLoc(String line) {
        return decode(line);
    }

    public void setSpliceJunctions(boolean b) {
        spliceJunctions = true;

    }

    public void setGffTags(boolean gffTags) {
        this.gffTags = gffTags;
    }

    public boolean isGffTags() {
        return this.gffTags;
    }
}
