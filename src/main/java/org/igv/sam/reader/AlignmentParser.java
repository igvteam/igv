
package org.igv.sam.reader;

import org.igv.sam.Alignment;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.IOException;

/**
 * @author jrobinso
 */
public interface AlignmentParser {

    Alignment readNextRecord(AsciiLineReader reader) throws IOException;

}
