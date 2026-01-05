
package org.igv.ui;

import org.igv.ui.util.CheckList;

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
