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

import net.sf.picard.sam.BuildBamIndex;
import net.sf.samtools.SAMFileReader;

import javax.swing.*;
import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Dec 6, 2009
 * Time: 7:42:01 PM
 * To change this template use File | Settings | File Templates.
 */
public class SamIndexer extends AlignmentIndexer {

    final static int FLAG_COL = 1;
    final static int READ_UNMAPPED_FLAG = 0x4;

    public SamIndexer(File samFile, JProgressBar progressBar, SamIndexCreatorDialog.IndexWorker worker) {
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

    //TODO MOVE THIS SOMEWHERE ELSE, LIKELY BamIndexer class
    //TODO Or maybe remove completely
    public static void createBAMIndex(File inputBam, File outputBam){
        SAMFileReader reader = new SAMFileReader(inputBam);
        BuildBamIndex.createIndex(reader, outputBam);
    }
}
