/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

package org.broad.igv.lists;

import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 * @date Sep 26, 2010
 */
public class GeneList {

    private String group;
    private boolean editable = true;
    private String name;
    private String description;
    private List<String> loci;

    public GeneList(String name, String description, String group, List<String> loci) {
        this.group = group;
        this.description = description;
        this.name = name;
        this.loci = loci;
    }

    public GeneList(String name, List<String> loci) {
        this.group = "My lists";
        this.name = name;
        this.loci = loci;
    }

    public String getName() {
        return name;
    }

    public List<String> getLoci() {
        return loci;
    }

    public int size() {
        return loci == null ? 0 : loci.size();
    }

    public void add(String gene) {
        // List might be immutable (Arrays.ArrayList)
        loci = new ArrayList(loci);
        loci.add(gene);
    }

    public GeneList copy() {
        return new GeneList(name + " copy", loci);
    }

    public void setLoci(List<String> strings) {
        loci = strings;
    }

    public void setName(String name) {
        this.name = name;
    }

    public boolean isEditable() {
        return editable;
    }

    public String getGroup() {
        
        return group;
    }

    public void setEditable(boolean editable) {
        this.editable = editable;
    }

    public void setGroup(String group) {
        this.group = group;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }
}
