package org.igv.util;

import org.igv.track.AttributeManager;
import org.igv.track.Track;
import org.igv.util.collections.CollUtils;

/**
 * @author eflakes
 */
public class FilterElement {

    public enum Operator implements CollUtils.Valued {

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

    public enum BooleanOperator {

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


    private String attributeKey;

    private Operator comparisonOperator = Operator.EQUAL;

    private String expectedValue;

    private BooleanOperator booleanOperator = BooleanOperator.AND;

    public FilterElement(String attributeKey, Operator comparisonOperator,
                         String expectedValue) {

        if (attributeKey != null) {
            this.attributeKey = attributeKey.trim();
        }

        if (comparisonOperator != null) {
            this.comparisonOperator = comparisonOperator;
        }

        if (expectedValue != null) {
            this.expectedValue = new String(expectedValue).toUpperCase().trim();
        }
    }

    @Override
        public String toString() {
            return "FilterElement{" +
                    "attributeKey='" + attributeKey + '\'' +
                    ", comparisonOperator=" + comparisonOperator +
                    ", expectedValue='" + expectedValue +
                    '}';
        }

    public boolean test(String attributeValue) {

        boolean result = false;

        // Changed by JTR.  Treat nulls as empty strings, so " == null" will
        // match only other nulls,   != null will match all non-nulls.
        if (expectedValue == null) {
            expectedValue = "";
        }

        if (attributeValue == null) {
            attributeValue = "";
        }

        // convert case
        attributeValue = attributeValue.toUpperCase().trim();

        boolean isComparableNumeric = isNumber(attributeValue);
        boolean isExpectedValueNumeric = isNumber(expectedValue);

        int comparison = 0;
        if (isComparableNumeric && isExpectedValueNumeric) {

            // Compare as numbers
            Double number1 = Double.parseDouble(attributeValue);
            Double number2 = Double.parseDouble(expectedValue);
            comparison = number1.compareTo(number2);
        } else if (isExpectedValueNumeric && attributeValue.equals("")) {
            // Null (blank) arguments are treated as < all other numbers when
            // comparing numercally
            Double number1 = Double.MIN_VALUE;
            Double number2 = Double.parseDouble(expectedValue);
            comparison = number1.compareTo(number2);

        } else {

            // Compare as strings
            comparison = attributeValue.compareTo(expectedValue);
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
            result = (attributeValue.indexOf(expectedValue) != -1);
        } else if (comparisonOperator.equals(Operator.DOES_NOT_CONTAIN)) {
            result = (attributeValue.indexOf(expectedValue) == -1);
        } else if (comparisonOperator.equals(Operator.STARTS_WITH)) {
            result = (attributeValue.startsWith(expectedValue));
        }

        return result;
    }

    public Operator getComparisonOperator() {
        return comparisonOperator;
    }

    public String getAttributeKey() {
        return attributeKey;
    }

    public String getValue() {
        return expectedValue;
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


    public boolean evaluate(Track track) {

        String attributeKey = getAttributeKey();
        String attribute = track.getAttributeValue(attributeKey);
        return test(attribute);
    }

    public boolean evaluateSample(String sample) {

        String attributeKey = getAttributeKey();
        String attribute = AttributeManager.getInstance().getAttribute(sample, attributeKey);

        return test(attribute);
    }
}
