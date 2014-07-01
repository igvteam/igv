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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import net.sf.samtools.util.CloseableIterator;

public class EmptyAlignmentIterator implements CloseableIterator<PicardAlignment> {

    static EmptyAlignmentIterator instance = new EmptyAlignmentIterator();

    public static EmptyAlignmentIterator getInstance() {
        return instance;
    }

    private EmptyAlignmentIterator() {
    }

    public void close() {
        // ignore
    }

    public boolean hasNext() {
        return false;
    }

    public PicardAlignment next() {
        return null;
    }

    public void remove() {
        // ignore
    }
}
