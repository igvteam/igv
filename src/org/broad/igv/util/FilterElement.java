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
package org.broad.igv.util;

import org.broad.igv.track.Track;
import org.broad.igv.util.collections.CollUtils;

/**
 * @author eflakes
 */
abstract public class FilterElement {

    public static enum Operator implements CollUtils.Valued {

        EQUAL("is equal to"),
        NOT_EQUAL("is not equal to"),
        GREATER_THAN("is greater than"),
        LESS_THAN("is less than"),
        GREATER_THAN_OR_EQUAL("is greater than or equal to"),
        LESS_THAN_OR_EQUAL("is less than or equal to"),
        STARTS_WITH("starts with"),
        CONTAINS("contains"),
        DOES_NOT_CONTAIN("does not contain");

        String value;

        Operator(String value) {
            this.value = value;
        }

        public String getValue() {
            return value;
        }
    }

    public static enum BooleanOperator {

        AND("AND"),
        OR("OR");

        String value;

        BooleanOperator(String value) {
            this.value = value;
        }

        public String getValue() {
            return value;
        }
    }

    private Filter filter;

    private String selectedItem;

    private Operator comparisonOperator = Operator.EQUAL;

    private String expectedValue;

    private BooleanOperator booleanOperator = BooleanOperator.AND;

    public FilterElement(Filter filter, String item, Operator comparisonOperator,
                         String expectedValue, BooleanOperator booleanOperator) {

        this.filter = filter;

        if (item != null) {
            this.selectedItem = item.trim();
        }

        if (comparisonOperator != null) {
            this.comparisonOperator = comparisonOperator;
        }

        if (expectedValue != null) {
            this.expectedValue = new String(expectedValue).toUpperCase().trim();
        }

        if (booleanOperator != null) {
            this.booleanOperator = booleanOperator;
        }

    }

    public boolean test(String comparableItem, Boolean previousResult) {

        boolean result = false;

        // Changed by JTR.  Treat nulls as empty strings, so " == null" will
        // match only other nulls,   != null will match all non-nulls.
        if (expectedValue == null) {
            expectedValue = "";
        }

        if (comparableItem == null) {
            comparableItem = "";
        }

        // convert case
        comparableItem = comparableItem.toUpperCase().trim();

        boolean isComparableNumeric = isNumber(comparableItem);
        boolean isExpectedValueNumeric = isNumber(expectedValue);

        int comparison = 0;
        if (isComparableNumeric && isExpectedValueNumeric) {

            // Compare as numbers
            Double number1 = Double.parseDouble(comparableItem);
            Double number2 = Double.parseDouble(expectedValue);
            comparison = number1.compareTo(number2);
        } else if (isExpectedValueNumeric && comparableItem.equals("")) {
            // Null (blank) arguments are treated as < all other numbers when
            // comparing numercally
            Double number1 = Double.MIN_VALUE;
            Double number2 = Double.parseDouble(expectedValue);
            comparison = number1.compareTo(number2);

        } else {

            // Compare as strings
            comparison = (comparableItem).compareTo(expectedValue);
        }

        if (comparisonOperator.equals(Operator.EQUAL)) {
            result = (comparison == 0);
        } else if (comparisonOperator.equals(Operator.NOT_EQUAL)) {
            result = (comparison != 0);
        } else if (comparisonOperator.equals(Operator.GREATER_THAN)) {
            result = (comparison > 0);
        } else if (comparisonOperator.equals(Operator.LESS_THAN)) {
            result = (comparison < 0);
        } else if (comparisonOperator.equals(Operator.GREATER_THAN_OR_EQUAL)) {
            result = (comparison >= 0);
        } else if (comparisonOperator.equals(Operator.LESS_THAN_OR_EQUAL)) {
            result = (comparison <= 0);
        } else if (comparisonOperator.equals(Operator.CONTAINS)) {
            result = (comparableItem.indexOf(expectedValue) != -1);
        } else if (comparisonOperator.equals(Operator.DOES_NOT_CONTAIN)) {
            result = (comparableItem.indexOf(expectedValue) == -1);
        } else if (comparisonOperator.equals(Operator.STARTS_WITH)) {
            result = (comparableItem.startsWith(expectedValue));
        }

        // If we have previous result we need to test against them
        if (previousResult != null) {

            // And/Or new result to previous result
            if (booleanOperator.equals(BooleanOperator.OR)) {
                result = (result || previousResult);
            } else {
                result = (result && previousResult);
            }
        }

        return result;
    }

    public void setBooleanOperator(BooleanOperator booleanOperator) {
        this.booleanOperator = booleanOperator;
    }

    public void setComparisonOperator(Operator comparisonOperator) {
        this.comparisonOperator = comparisonOperator;
    }

    public void setSelectedItem(String item) {
        this.selectedItem = item;
    }

    public void setExpectedValue(String expectedValue) {
        if (expectedValue != null) {
            this.expectedValue = new String(expectedValue).toUpperCase();
        }
    }

    public BooleanOperator getBooleanOperator() {
        return booleanOperator;
    }

    public Operator getComparisonOperator() {
        return comparisonOperator;
    }

    public String getSelectedItem() {
        return selectedItem;
    }

    public String getValue() {
        return expectedValue;
    }

    public Filter getFilter() {
        return filter;
    }

    public boolean isNumber(String string) {
        if (string == null) {
            return false;
        }

        string = string.trim();
        if (string.equals("")) {
            return false;
        }

        char charcters[] = string.toCharArray();
        boolean withDecimal = false, isNegative = false;
        int len = charcters.length;

        for (int i = 0; i < len; i++) {
            if (!Character.isDigit(charcters[i])) {

                if (charcters[i] == '.') {

                    if (withDecimal) {
                        return false;
                    }
                    withDecimal = true;
                } else if (charcters[i] == '-') {

                    if (isNegative) {
                        return false;
                    }

                    if (i == 0) {
                        isNegative = true;
                    } else {
                        return false;
                    }

                } else {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @param previousResult Result of the previous FilterElement that is
     *                       being chained to this one (null if no previous FilterElement);
     * @return
     */
    abstract public boolean evaluate(Track track, Boolean previousResult);
}
