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

package org.broad.igv.data.cufflinks;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Range;

/**
 * @author jrobinso
 *         Date: 3/8/13
 *         Time: 9:32 PM
 */
abstract public class CufflinksValue extends Range implements LocusScore {
    String gene;

    public CufflinksValue(String chr, int start, int end, String gene) {
        super(chr, start, end);
        this.gene = gene;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public String getGene() {
        return gene;
    }

}
