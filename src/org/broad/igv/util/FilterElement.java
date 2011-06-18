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
package org.broad.igv.util;

import org.broad.igv.track.Track;

/**
 * @author eflakes
 */
abstract public class FilterElement {

    public static enum Operator {

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

        static public Operator findEnum(String value) {

            if (value == null) {
                return null;
            }

            if (value.equals(EQUAL.getValue())) {
                return EQUAL;
            } else if (value.equals(NOT_EQUAL.getValue())) {
                return NOT_EQUAL;
            } else if (value.equals(GREATER_THAN.getValue())) {
                return GREATER_THAN;
            } else if (value.equals(LESS_THAN.getValue())) {
                return LESS_THAN;
            } else if (value.equals(GREATER_THAN_OR_EQUAL.getValue())) {
                return GREATER_THAN_OR_EQUAL;
            } else if (value.equals(LESS_THAN_OR_EQUAL.getValue())) {
                return LESS_THAN_OR_EQUAL;
            } else if (value.equals(STARTS_WITH.getValue())) {
                return STARTS_WITH;
            } else if (value.equals(CONTAINS.getValue())) {
                return CONTAINS;
            } else if (value.equals(DOES_NOT_CONTAIN.getValue())) {
                return CONTAINS;
            }
            return null;
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

        static public BooleanOperator findEnum(String value) {

            if (value == null) {
                return null;
            }

            if (value.equals(AND.getValue())) {
                return AND;
            } else if (value.equals(OR.getValue())) {
                return OR;
            }
            return null;
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
