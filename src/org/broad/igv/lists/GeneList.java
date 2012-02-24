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

package org.broad.igv.lists;

import org.broad.igv.feature.Locus;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
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
    private boolean showName = true;
    private static Comparator<String> POSITION_COMPARATOR;

    public GeneList(String name, String description, String group, List<String> loci) {
        this.group = group;
        this.description = description;
        this.name = name;
        //We do this to guarantee that certain operations will be supported
        this.loci = loci;
    }

    public GeneList(String name, List<String> loci) {
        this(name, null, GeneListManager.USER_GROUP, loci);
    }

    public GeneList(String name, List<String> loci, boolean showName) {
        this(name, null, GeneListManager.USER_GROUP, loci);
        this.showName = showName;
    }

    public GeneList() {
        this.group = GeneListManager.USER_GROUP;
    }

    public String getName() {
        return name;
    }


    public String getDisplayName() {
        return showName ? name : "";
    }

    public List<String> getLoci() {
        return loci;
    }

    public int size() {
        return loci == null ? 0 : loci.size();
    }

    public void add(String gene) {
        if (loci == null) {
            loci = new ArrayList<String>(1);
        }
        try {
            //Can't guarantee that list will support this operation
            loci.add(gene);
        } catch (Exception e) {
            loci = new ArrayList<String>(loci);
            loci.add(gene);
        }
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


    /**
     * Sort loci by "position".  This only sorts loci of the form chr1:100-200.
     */
    public static void sortByPosition(List<String> loci) {
        if (POSITION_COMPARATOR == null) initComparator();

        Collections.sort(loci, POSITION_COMPARATOR);
    }


    private static synchronized void initComparator() {
        POSITION_COMPARATOR = new Comparator<String>() {
            public int compare(String s1, String s2) {
                Locus l1 = new Locus(s1);
                Locus l2 = new Locus(s2);
                if (!l1.isValid() && !l2.isValid()) {
                    return 0;
                } else if (!l1.isValid()) {
                    return -1;
                } else if (!l2.isValid()) {
                    return 1;
                } else if (!l1.getChr().equals(l2.getChr())) {
                    return l1.getChr().compareTo(l2.getChr());
                } else {
                    return l1.getStart() - l2.getStart();
                }
            }

        };
    }

}
