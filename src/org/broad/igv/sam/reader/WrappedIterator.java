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

package org.broad.igv.sam.reader;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.sam.PicardAlignment;

/**
 * @author jrobinso
 * @since Sep 22, 2009
 */
public class WrappedIterator implements CloseableIterator<PicardAlignment> {

    CloseableIterator<SAMRecord> iter;

    public WrappedIterator(CloseableIterator<SAMRecord> iter) {
        this.iter = iter;
    }

    public void close() {
        iter.close();
    }

    public boolean hasNext() {
        return iter.hasNext();
    }

    public PicardAlignment next() {
        return new PicardAlignment(iter.next());
    }

    public void remove() {
        iter.remove();
    }
}
