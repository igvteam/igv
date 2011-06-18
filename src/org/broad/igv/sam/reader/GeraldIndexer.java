/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.sam.reader;

import javax.swing.*;
import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Dec 6, 2009
 * Time: 7:42:29 PM
 * To change this template use File | Settings | File Templates.
 */ /*
*         String readChr = mapChromosome(fields[CHROMOSOME_COLUMN]);
int alignmentStart = Integer.parseInt(fields[ALIGNMENT_START_COLUMN].trim()) - 1;
*/
public class GeraldIndexer extends AlignmentIndexer {

    public GeraldIndexer(File samFile, JProgressBar progressBar, SamIndexCreatorDialog.IndexWorker worker) {
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
