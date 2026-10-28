
package org.igv.alignment.reader;

import org.igv.alignment.Alignment;
import htsjdk.tribble.readers.Utf8LineReader;

import java.io.IOException;

/**
 * @author jrobinso
 */
public interface AlignmentParser {

    Alignment readNextRecord(Utf8LineReader reader) throws IOException;

}
