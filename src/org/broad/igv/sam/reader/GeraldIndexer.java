/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
public class GeraldIndexer extends AlignmentIndexer {

    public GeraldIndexer(File samFile, JProgressBar progressBar, IndexCreatorDialog.SamIndexWorker worker) {
        super(samFile, progressBar, worker);
    }

    int getAlignmentStart(String[] fields) throws NumberFormatException {
        // Get alignmetn start and verify file is sorted.
        return Integer.parseInt(fields[GeraldParser.ALIGNMENT_START_COLUMN].trim()) - 1;
    }

    int getAlignmentLength(String[] fields) throws NumberFormatException {
        return fields[GeraldParser.READ_COLUMN].length();
    }

    String getChromosome(String[] fields) {
        return fields[GeraldParser.CHROMOSOME_COLUMN];
    }

    @Override
    boolean isMapped(String[] fields) {
        return true;
    }
}
