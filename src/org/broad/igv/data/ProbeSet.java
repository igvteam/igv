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

/*
 * ProbeSet.java
 *
 * Created on October 31, 2007, 5:18 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package org.broad.igv.data;

import java.util.*;

/**
 * @author jrobinso
 */
public class ProbeSet {

    private Map<String, List<ExpressionProbe>> probeLists = new HashMap();
    private Map<String, ExpressionProbe> probes = new HashMap();

    public Set<String> getChromosomes() {
        return probeLists.keySet();
    }


    public List<ExpressionProbe> getProbes(String chr) {
        return probeLists.get(chr);
    }

    public ExpressionProbe getProbe(String probeId) {
        return probes.get(probeId);
    }

    public void add(ExpressionProbe probe) {
        probes.put(probe.getName(), probe);
        List<ExpressionProbe> pList = probeLists.get(probe.getChr());
        if (pList == null) {
            pList = new ArrayList(1000);
            probeLists.put(probe.getChr(), pList);
        }
        pList.add(probe);
    }

    public void sortProbeLists() {
        for (List<ExpressionProbe> list : probeLists.values()) {
            Collections.sort(list);
        }
    }


}
