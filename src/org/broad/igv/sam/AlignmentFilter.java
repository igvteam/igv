package org.broad.igv.sam;

import net.sf.samtools.SAMRecord;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.broad.igv.util.FilterElement.BooleanOperator;
import org.broad.igv.util.FilterElement.Operator;

public class AlignmentFilter {

	public boolean isExclude() {
		return exclude;
	}

	public void setExclude(boolean exclude) {
		this.exclude = exclude;
	}

	public String getExpression() {
		return expression;
	}

	public void setExpression(String expression) {
		this.expression = expression;
	}

	public void setTag(String tag) {
		this.tag = tag;
	}

	public static enum Operator {

		EQUAL("is equal to"), NOT_EQUAL("is not equal to"), GREATER_THAN(
				"is greater than"), LESS_THAN("is less than"), GREATER_THAN_OR_EQUAL(
				"is greater than or equal to"), LESS_THAN_OR_EQUAL(
				"is less than or equal to"), STARTS_WITH("starts with"), CONTAINS(
				"contains"), REGULAR_EXPRESSION("is regular expression"), CONTAINS_WILDCARDS(
				"contains wild cards"), IS_MISSING("is missing"), DOES_NOT_CONTAIN(
				"does not contain");

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
			} else if (value.equals(REGULAR_EXPRESSION.getValue())) {
				return REGULAR_EXPRESSION;
			} else if (value.equals(CONTAINS_WILDCARDS.getValue())) {
				return CONTAINS_WILDCARDS;
			} else if (value.equals(IS_MISSING.getValue())) {
				return IS_MISSING;
			} else if (value.equals(DOES_NOT_CONTAIN.getValue())) {
				return CONTAINS;
			}
			return null;
		}
	}

	private Operator comparisonOperator = Operator.EQUAL;

	public static enum BooleanOperator {

		AND("AND"), OR("OR");

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

	private boolean exclude;
	private String tag;
	private String expression;
	private String Operation;

	AlignmentFilter() {
		exclude = true;
		tag = "NM";
		expression = "*T*";
		Operation = "does not contain";
	}

	public String getTag() {
		return tag;
	}

	//false = keep for display
	//true  = filter out (remove from display)
	public boolean filter(String attr) {
		String regExprToUse = null;
		Pattern regExpr;

		boolean match = false;

		if (Operation.equals("is equal to")) {
			if (attr.equals(expression)) {
				match = true;
			}
		} else if (Operation.equals("is not equal to")) {
			if (attr.equals(expression)) {
				match = true;
			}
		} else if (Operation.equals("starts with")) {
			if (attr.startsWith(expression)) {
				match = true;
			}
		} else if (Operation.equals("contains")) {
			if (attr.contains(expression)) {
				match = true;
			}
		} else if (Operation.equals("is regular expression")) {
			regExprToUse = expression;
			regExpr = Pattern.compile(regExprToUse);
			Matcher matcher = regExpr.matcher(attr);
			match = matcher.matches();
		} else if (Operation.equals("does not contain")) {
			if (attr.contains(expression)) {
				match = true;
			}
		} else if (Operation.equals("is missing")) {
			//if the flag is missing (that is what "is missing" is referring to, we won't be here
			return false;
		}

		//if (hasWildCards) {
			// regExprToUse = WildcardMatcher.wildcardToRegex(expression);
		//} else 
	

		return match;
	}
	public boolean filter(Double attr) {
		if (Operation.equals("is greater than")) {
			return Double.parseDouble(expression) > attr;
		} else if (Operation.equals("is less than")) {
			return Double.parseDouble(expression) < attr;
		} else if (Operation.equals("is greater than or equal to")) {
			return Double.parseDouble(expression) >= attr;
		} else if (Operation.equals("is less than or equal to")) {
			return Double.parseDouble(expression) <= attr;
		} else{
			return filter(attr.toString());
		}
	}
	public boolean filter(Object attribute) {
		if (attribute instanceof Integer) {
			return filter(new Double((Integer)attribute));
		}else{
			return filter(attribute.toString());
		}
	}

	public void setOperation(String op) {
		this.Operation = op;
	}

	public String getOperation() {
		return this.Operation;
	}

	// public boolean filter(Object attr) {
	// return false;
	// }

}
