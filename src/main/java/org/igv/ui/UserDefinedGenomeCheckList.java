/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.igv.ui;

import org.igv.ui.util.CheckList;

import java.util.HashSet;

/**
 * @author eflakes
 */
public class UserDefinedGenomeCheckList extends CheckList {

    public UserDefinedGenomeCheckList(boolean defaultState) {
        super(defaultState, "<html>Select one or more imported genomes to remove.</html>");
    }

    public HashSet<String> getSelectedGenomes() {
        return getSelectedItems();
    }

    public HashSet<String> getUnselectedGenomes() {
        return getUnselectedItems();
    }
}
