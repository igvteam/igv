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

import org.broad.igv.ui.util.IndexCreatorDialog;

import javax.swing.*;
import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * @author jrobinso
 * @since Dec 6, 2009
 */
public class DotAlignedIndexer extends AlignmentIndexer {

    int baseOffset = 1;

    public DotAlignedIndexer(File samFile, JProgressBar progressBar, IndexCreatorDialog.SamIndexWorker worker) {
        super(samFile, progressBar, worker);
        if (samFile.getName().endsWith(".bedz") || samFile.getName().endsWith(".bed")) {
            baseOffset = 0;
        }
    }

    int getAlignmentStart(String[] fields) throws NumberFormatException {
        int position = Integer.parseInt(fields[1]) - baseOffset;
        return position;
    }

    int getAlignmentLength(String[] fields) throws NumberFormatException {
        return Integer.parseInt(fields[2]) - Integer.parseInt(fields[1]) + 1;
    }

    String getChromosome(String[] fields) {
        String chr = fields[0];
        return chr;
    }

    @Override
    boolean isMapped(String[] fields) {
        return true;
    }
}
