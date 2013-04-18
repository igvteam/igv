/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.bbfile;

import net.sf.samtools.seekablestream.SeekableStream;

/**
 * @author jrobinso
 * @date Jun 22, 2011
 */
public class RPTreeNodeProxy  {

    public SeekableStream fis;
    public long fileOffset;
    public boolean isLowToHigh;

    // For debugging
    int chromId;

    public RPTreeNodeProxy(SeekableStream fis, long fileOffset, boolean lowToHigh, int chromId) {
        this.fis = fis;
        this.fileOffset = fileOffset;
        isLowToHigh = lowToHigh;
        this.chromId = chromId;
    }

}
