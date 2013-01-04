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

package org.broad.igv.session;

import org.apache.log4j.Logger;
import org.broad.igv.lists.GeneListManager;
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
                IGV.getInstance().setGeneList(GeneListManager.getInstance().getGeneList(listName), false);

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
