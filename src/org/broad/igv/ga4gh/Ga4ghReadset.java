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

package org.broad.igv.ga4gh;

/**
 * Minimal representation of a readset for prototype purposes.
 */
public class Ga4ghReadset {

    String id;
    String name;
    String genomeId;

    public Ga4ghReadset(String id, String name, String genomeId) {
        this.id = id;
        this.name = name;
        this.genomeId = genomeId;
    }

    public String getId() {
        return id;
    }

    public String getName() {
        return name;
    }

    public String getGenomeId() {
        return genomeId;
    }

    public String toString() {
        return name;
    }
}
