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

package org.broad.igv.ui;

import javax.swing.*;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

class JList7<T> extends JList {

    public List<T> getSelectedValuesList() {
        ArrayList<T> typedSelected = new ArrayList<T>(getSelectedValues().length);
        for (Object selected : getSelectedValues()) {
            try {
                typedSelected.add((T) selected);
            } catch (ClassCastException exc) {
                //pass
            }
        }
        return typedSelected;
    }

}