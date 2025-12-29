package org.broad.igv.sam.reader;

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
