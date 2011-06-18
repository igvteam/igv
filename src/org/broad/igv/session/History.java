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

package org.broad.igv.session;

import org.apache.log4j.Logger;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.action.SearchCommand;
import org.broad.igv.ui.panel.FrameManager;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * @author jrobinso
 * @date Sep 7, 2010
 */
public class History {

    private static Logger log = Logger.getLogger(History.class);

    int maxEntries = 100;
    int currPos = 0;

    LinkedList<Entry> activeStack;
    List<Entry> allHistory;

    public History(int maxEntries) {
        this.maxEntries = maxEntries;
        activeStack = new LinkedList();
        allHistory = new ArrayList();
    }


    public void push(String s, int zoom) {

        if (s == null || s.length() == 0) {
            return;
        }
        // If in gene list mode disable history,  with the exception of "List" items
        if (FrameManager.isGeneListMode() && !s.startsWith("List")) {
            return;
        }

       if(activeStack.size() > 0 && s.equals(activeStack.peek().locus)) {
            return;
        }


        log.debug("History: " + s);
        allHistory.add(new Entry(s, zoom));

        while (currPos > 0) {
            activeStack.removeFirst();
            currPos--;
        }
        activeStack.addFirst(new Entry(s, zoom));
    }


    public void back() {
        if (canGoBack()) {
            currPos++;
            Entry item = activeStack.get(currPos);
            processItem(item);
        }
    }

    public void forward() {

        if (canGoForward()) {
            currPos--;
            Entry item = activeStack.get(currPos);
            processItem(item);
        }
    }

    public boolean canGoBack() {
        return activeStack.size() > 0 && currPos < (activeStack.size() - 1);
    }

    public boolean canGoForward() {
        return activeStack.size() > 0 && currPos > 0;
    }

    public void processItem(Entry entry) {

        if (entry != null) {
            String locus = entry.getLocus();
            if (locus.startsWith("List: ")) {
                String listName = locus.substring(6);
                IGV.getInstance().setGeneList(listName, false);

            } else {
                if (FrameManager.isGeneListMode()) {
                    IGV.getInstance().setGeneList(null, false);
                }
                (new SearchCommand(FrameManager.getDefaultFrame(), locus, false)).execute();
                FrameManager.getDefaultFrame().setZoom(entry.getZoom());
            }
        }
    }


    public Entry peekBack() {
        if (activeStack.size() == 0 || (currPos + 1) >= activeStack.size()) {
            return null;
        }
        return activeStack.get(currPos + 1);
    }

    public Entry peekForward() {
        if (activeStack.size() == 0 || (currPos - 1) < 0) {
            return null;
        }
        return activeStack.get(currPos - 1);
    }

    public void clear() {
        activeStack.clear();
        allHistory.clear();
        currPos = 0;
    }

    public void printStack() {
        System.out.println("curr pos=" + currPos);
        for (Entry s : activeStack) {
            System.out.println(s.getLocus());
        }
        System.out.println();
    }

    public List<Entry> getAllHistory() {
        return allHistory;
    }

    public static class Entry {
        private String locus;
        private int zoom;

        public Entry(String s, int zoom) {
            this.locus = s;
            this.zoom = zoom;
        }

        public String getLocus() {
            return locus;
        }

        public int getZoom() {
            return zoom;
        }

        public String toString() {
            return locus;
        }
    }

}
