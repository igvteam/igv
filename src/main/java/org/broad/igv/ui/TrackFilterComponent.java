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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.ui;

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

    public FilterElement createFilterElement(Filter filter, String selectedItem,
                                                  Operator comparisonOperator, String value, BooleanOperator booleanOperator) {

        FilterElement filterElement =
                new FilterElement(
                        filter,
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
