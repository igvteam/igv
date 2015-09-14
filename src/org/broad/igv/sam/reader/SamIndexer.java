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

package org.broad.igv.sam.reader;

import picard.sam.BuildBamIndex;
import htsjdk.samtools.SAMFileReader;
import org.broad.igv.ui.util.IndexCreatorDialog;

import javax.swing.*;
import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * @author jrobinso
 * @since: Dec 6, 2009
 */
public class SamIndexer extends AlignmentIndexer {

    final static int FLAG_COL = 1;
    final static int READ_UNMAPPED_FLAG = 0x4;

    public SamIndexer(File samFile, JProgressBar progressBar, IndexCreatorDialog.SamIndexWorker worker) {
        super(samFile, progressBar, worker);
    }

    int getAlignmentStart(String[] fields) throws NumberFormatException {
        // Get alignmetn start and verify file is sorted.
        int alignmentStart = Integer.parseInt(fields[3].trim()) - 1;
        return alignmentStart;
    }

    int getAlignmentLength(String[] fields) throws NumberFormatException {
        String cigarString = fields[5];
        return SamUtils.getPaddedReferenceLength(cigarString);
    }

    String getChromosome(String[] fields) {
        return fields[2];
    }

    @Override
    boolean isMapped(String[] fields) {
        int flags = Integer.parseInt(fields[FLAG_COL]);
        return (flags & READ_UNMAPPED_FLAG) == 0;
    }

}
