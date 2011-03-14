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

package org.broad.igv.session;

import org.apache.log4j.Logger;
import org.broad.igv.ui.IGVMainFrame;
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

    LinkedList<String> activeStack;
    List<String> allHistory;

    public History(int maxEntries) {
        this.maxEntries = maxEntries;
        activeStack = new LinkedList();
        allHistory = new ArrayList();
    }


    public void push(String s) {

        // If in gene list mode disable history,  with the exception of "List" items

        if (FrameManager.isGeneListMode() && !s.startsWith("List")) {
            return;
        }

        log.debug("History: " + s);
        allHistory.add(s);

        while (currPos > 0) {
            activeStack.removeFirst();
            currPos--;
        }
        activeStack.addFirst(s);


    }


    public void back() {
        if (activeStack.size() == 0 || currPos >= (activeStack.size() - 1)) {
            return;
        }
        currPos++;
        String item = activeStack.get(currPos);
        processItem(item);
    }

    public void forward() {

        if (activeStack.size() == 0 || currPos == 0) {
            return;
        }
        currPos--;
        printStack();
        String item = activeStack.get(currPos);
        processItem(item);
    }

    public void processItem(String item) {
        if (item != null) {
            if (item.startsWith("List: ")) {
                String listName = item.substring(6);
                IGVMainFrame.getInstance().setGeneList(listName, false);

            } else {
                if (FrameManager.isGeneListMode()) {
                    IGVMainFrame.getInstance().setGeneList(null, false);
                }
                (new SearchCommand(FrameManager.getDefaultFrame(), item, false)).execute();
                IGVMainFrame.getInstance().refreshCommandBar();
            }
        }
    }


    public String peekBack() {
        if (activeStack.size() == 0 || (currPos + 1) >= activeStack.size()) {
            return null;
        }
        return activeStack.get(currPos + 1);
    }

    public String peekForward() {
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
        for (String s : activeStack) {
            System.out.println(s);
        }
        System.out.println();
    }

    public List<String> getAllHistory() {
        return allHistory;
    }

}
