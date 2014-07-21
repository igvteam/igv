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

package org.broad.igv.feature.tribble;

import htsjdk.tribble.Feature;

/**
 * A minimal representation of a Feature.
 */
public class Locus implements Feature {

    String chr;
    int start;
    int end;

    /**
     *
     * @param chr
     * @param start
     * @param end
     */
    public Locus(String chr, int start, int end) {
        this.chr = chr;
        this.start = start;
        this.end = end;
    }

    public String getChr() {
        return chr;
    }

    public int getStart() {
        return start;  
    }

    public int getEnd() {
        return end;
    }
}
