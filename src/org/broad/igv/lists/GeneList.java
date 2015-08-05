/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.lists;

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
        init(name, description, group, loci);
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

    private void init(String name, String description, String group, List<String> loci) {
        this.group = group;
        this.description = description;
        this.name = name;
        //We do this to guarantee that certain operations will be supported
        this.loci = loci;
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
//    public static void sortByPosition(List<String> loci) {
//        if (POSITION_COMPARATOR == null) initComparator();
//
//        Collections.sort(loci, POSITION_COMPARATOR);
//    }
//
//    private static synchronized void initComparator() {
//        POSITION_COMPARATOR = new Comparator<String>() {
//            public int compare(String s1, String s2) {
//                Locus l1 = Locus.fromString(s1);
//                Locus l2 = Locus.fromString(s2);
//                boolean l1Valid = l1 != null;
//                boolean l2Valid = l2 != null;
//                if (l1Valid && l2Valid) {
//                    return 0;
//                } else if (!l1Valid) {
//                    return -1;
//                } else if (!l2Valid) {
//                    return 1;
//                } else if (!l1.getChr().equals(l2.getChr())) {
//                    return l1.getChr().compareTo(l2.getChr());
//                } else {
//                    return l1.getStart() - l2.getStart();
//                }
//            }
//
//        };
//    }

    public void sort(Comparator<String> comparator) {
        Collections.sort(this.loci, comparator);
    }
}
