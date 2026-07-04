
package org.igv.alignment.reader;

import org.igv.alignment.Alignment;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.IOException;

/**
 * @author jrobinso
 */
public interface AlignmentParser {

    Alignment readNextRecord(AsciiLineReader reader) throws IOException;

}
