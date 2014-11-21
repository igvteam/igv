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

import java.util.List;

/**
 * Created by jrobinso on 9/3/14.
 */
public class Ga4ghDataset {

    String id;
    String name;
    String genomeId;
    List<Ga4ghReadset> readsets;

    public Ga4ghDataset(String id, String name, String genomeId) {
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

    public List<Ga4ghReadset> getReadsets() {
        return readsets;
    }

    public void setReadsets(List<Ga4ghReadset> readsets) {
        this.readsets = readsets;
    }
}
