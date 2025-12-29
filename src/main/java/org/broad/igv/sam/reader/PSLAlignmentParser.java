package org.broad.igv.sam.reader;


import org.broad.igv.feature.BasicFeature;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.FeatureWrappedAlignment;
import org.broad.igv.feature.tribble.PSLCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.IOException;

/**
 * @author jrobinso
 * @date Aug 5, 2010
 */

public class PSLAlignmentParser implements AlignmentParser {

    PSLCodec codec;

    public PSLAlignmentParser() {
        codec = new PSLCodec(); //genome);
    }

    public Alignment readNextRecord(AsciiLineReader reader) throws IOException {

        String nextLine;
        while ((nextLine = reader.readLine()) != null) {

            Feature f = codec.decode(nextLine);
            if(f != null && f instanceof BasicFeature) {
                return new FeatureWrappedAlignment((BasicFeature) f);
            }
        }
        return null;

    }
}
