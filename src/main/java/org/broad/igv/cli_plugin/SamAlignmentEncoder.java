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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.SAMTextWriter;
import org.broad.igv.sam.PicardAlignment;

import java.io.OutputStream;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Encode an Alignment into SAM format
 * @author jacob
 * @since 2012-Sep-27
 */
public class SamAlignmentEncoder implements FeatureEncoder<PicardAlignment> {

    private boolean headerSet = false;

    public Map<String, Object> encodeAll(OutputStream stream, Iterator<? extends PicardAlignment> alignments) {
        SAMFileWriterImpl writer = new SAMTextWriter(stream);
        while (alignments.hasNext()) {
            PicardAlignment samAl = alignments.next();
            if (!headerSet) {
                writer.setSortOrder(SAMFileHeader.SortOrder.unsorted, true);
                writer.setHeader(samAl.getRecord().getHeader());
                headerSet = true;
            }
            writer.addAlignment(samAl.getRecord());
        }
        writer.close();
        return null;
    }



    @Override
    public void setInputs(List<String> commands, Map<Argument, Object> argumentMap, Argument argument) {
        //pass
    }
}
