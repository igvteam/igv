/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.ui;

import org.broad.igv.ui.util.CheckList;

import java.util.HashSet;

/**
 * @author eflakes
 */
public class AttributeCheckList extends CheckList {

    public AttributeCheckList(boolean defaultState) {
        super(defaultState);
    }


    public HashSet<String> getUnselectedAttributes() {
        return getUnselectedItems();
    }
}
