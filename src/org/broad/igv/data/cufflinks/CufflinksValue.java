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

/**
 * @author jrobinso
 *         Date: 3/8/13
 *         Time: 9:32 PM
 */
abstract public class CufflinksValue implements LocusScore {
    String chr;
    int start;
    int end;

    String gene;

    public CufflinksValue(String gene, String chr, int start, int end) {
        this.gene = gene;
        this.chr = chr;
        this.start = start;
        this.end = end;
    }

    @Override
    public String getChr() {
        return chr;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getEnd() {
        return end;
    }

    @Override
    public void setStart(int start) {
        this.start = start;
    }

    @Override
    public void setEnd(int end) {
        this.end = end;
    }

    public String getGene() {
        return gene;
    }

}
