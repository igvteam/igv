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

package org.broad.igv.feature.tribble;

import org.broad.igv.renderer.SpliceJunctionRenderer;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.readers.LineReader;

import java.io.IOException;

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
    public Object readHeader(LineReader reader) {
        String line;
        try {
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#") || line.startsWith("track") ||
                        line.startsWith("browser")) {
                    readHeaderLine(line);
                } else {
                    break;
                }
            }
            return header;
        } catch (IOException e) {
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
            gffTags = tp.isGffTags();

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

}
