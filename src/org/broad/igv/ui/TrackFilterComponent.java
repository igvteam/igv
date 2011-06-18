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
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.ui;

import org.broad.igv.ui.util.FilterComponent;
import org.broad.igv.util.Filter;
import org.broad.igv.util.FilterElement;
import org.broad.igv.util.FilterElement.BooleanOperator;
import org.broad.igv.util.FilterElement.Operator;

import java.util.List;

/**
 * @author eflakes
 */
public class TrackFilterComponent extends FilterComponent {

    private boolean matchAll = true;

    public TrackFilterComponent(TrackFilterPane filterPane, String text, List<String> items,
                                FilterElement element) {

        super(filterPane, text, items, element);
    }

    public TrackFilterElement createFilterElement(Filter filter, String selectedItem,
                                                  Operator comparisonOperator, String value, BooleanOperator booleanOperator) {

        TrackFilterElement filterElement =
                new TrackFilterElement(
                        (TrackFilter) filter,
                        selectedItem,
                        comparisonOperator,
                        value,
                        booleanOperator);

        return filterElement;
    }

    public void setMatchAll(boolean value) {
        matchAll = value;
    }

    public boolean getMatchAll() {
        return matchAll;
    }

    /**
     * Save the UI content into a non-UI version of the FilterElement
     */
    public void save() {

        FilterElement filterElement = getFilterElement();

        // Boolean operator
        if (matchAll) {
            filterElement.setBooleanOperator(BooleanOperator.AND);
        } else {
            filterElement.setBooleanOperator(BooleanOperator.OR);
        }

        super.save();
    }

}
