/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

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
