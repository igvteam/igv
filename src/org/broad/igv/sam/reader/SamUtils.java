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

//~--- non-JDK imports --------------------------------------------------------

//~--- JDK imports ------------------------------------------------------------

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.FileUtils;

import java.io.File;
import java.io.IOException;

/**
 * @author jrobinso
 */
public class SamUtils {

    private static final byte ZERO_BYTE = "0".getBytes()[0];
    private static final byte NINE_BYTE = "9".getBytes()[0];

    public static FeatureIndex getIndexFor(String samPath) {

        String idxPath = samPath + ".sai";

        if (FileUtils.resourceExists(idxPath)) {
            return new FeatureIndex(idxPath);
        } else if (FileUtils.isRemote(idxPath)) {
            return null;
        } else {
            File idxFile = new File(idxPath);
            File samFile = new File(samPath);
            if (!idxFile.exists()) {
                idxFile = getUserIdxFile(samFile);
            }
            if (idxFile.exists() && idxFile.lastModified() > samFile.lastModified()) {
                return new FeatureIndex(idxFile);
            } else {
                return createIndexFor(samFile);
            }
        }
    }

    private static FeatureIndex createIndexFor(File samFile) {
        File newIdxFile = new File(samFile.getAbsolutePath() + ".sai");
        if (!FileUtils.canWriteTo(newIdxFile)) {
            newIdxFile = getUserIdxFile(samFile);
        }

        if (!Globals.isHeadless()) {
            SamIndexCreatorDialog dialog = new SamIndexCreatorDialog(IGV.getMainFrame(), true, samFile, newIdxFile);
            dialog.setLocationRelativeTo(IGV.getMainFrame());
            dialog.setVisible(true);
            return dialog.getIndex();
        } else {
            AlignmentIndexer indexer = AlignmentIndexer.getInstance(samFile, null, null);
            FeatureIndex index = null;
            try {
                System.out.println("Creating index");
                index = indexer.createSamIndex(newIdxFile);
            } catch (IOException e) {
                e.fillInStackTrace();
            }
            return index;
        }
    }

    private static File getUserIdxFile(File samFile) {
        File idxFile;
        File samDir = DirectoryManager.getSamDirectory();
        //Need the path information to distinguish like name indices in separate
        // directories.
        idxFile = new File(samDir, samFile.getName() + "_" + samFile.getParent().hashCode() + ".sai");
        return idxFile;
    }

    public static int getPaddedReferenceLength(String cigarString) {
        return decodeCigar(cigarString).getPaddedReferenceLength();

    }

    /**
     * Convert from String CIGAR representation to Cigar class representation.  Does not
     * do validation beyond the most basic CIGAR string well-formedness, i.e. each operator is
     * valid, and preceded by a decimal length.
     *
     * @param textCigar CIGAR in String form ala SAM text file.  "*" means empty CIGAR.
     * @throws RuntimeException if textCigar is invalid at the most basic level.
     */
    static Cigar decodeCigar(final String textCigar) {
        if (SAMRecord.NO_ALIGNMENT_CIGAR.equals(textCigar)) {
            return new Cigar();
        }
        final Cigar ret = new Cigar();
        final byte[] cigarBytes = StringUtil.stringToBytes(textCigar);
        for (int i = 0; i < cigarBytes.length; ++i) {
            if (!isDigit(cigarBytes[i])) {
                throw new IllegalArgumentException("Malformed CIGAR string: " + textCigar);
            }
            int length = (cigarBytes[i] - ZERO_BYTE);
            for (++i; isDigit(cigarBytes[i]); ++i) {
                length = (length * 10) + cigarBytes[i] - ZERO_BYTE;
            }
            final CigarOperator operator = CigarOperator.characterToEnum(cigarBytes[i]);
            ret.add(new CigarElement(length, operator));
        }
        return ret;
    }

    private static boolean isDigit(final byte c) {
        return c >= ZERO_BYTE && c <= NINE_BYTE;
    }

}
