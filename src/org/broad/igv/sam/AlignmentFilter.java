package org.broad.igv.sam;

import net.sf.samtools.SAMRecord;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class AlignmentFilter {

	public boolean isExclude() {
		return exclude;
	}

	public void setExclude(boolean exclude) {
		this.exclude = exclude;
	}

	public FilterType getFilterType() {
		return filterType;
	}

	public void setFilterType(FilterType filterType) {
		this.filterType = filterType;
	}

	public Double getMin() {
		return min;
	}

	public void setMin(Double min) {
		this.min = min;
	}

	public Double getMax() {
		return max;
	}

	public void setMax(Double max) {
		this.max = max;
	}

	public String getExpression() {
		return expression;
	}

	public void setExpression(String expression) {
		this.expression = expression;
	}

	public boolean isHasWildCards() {
		return hasWildCards;
	}

	public void setHasWildCards(boolean hasWildCards) {
		this.hasWildCards = hasWildCards;
	}

	public boolean isRegExpr() {
		return isRegExpr;
	}

	public void setRegExpr(boolean isRegExpr) {
		this.isRegExpr = isRegExpr;
	}

	public boolean isCaseSensitive() {
		return caseSensitive;
	}

	public void setCaseSensitive(boolean caseSensitive) {
		this.caseSensitive = caseSensitive;
	}

	public void setTag(String tag) {
		this.tag = tag;
	}

	public enum FilterType {
		// filter for a numeric range
		range,
		// filter for a string pattern
		pattern,
		// filter for missing values
		missing
	}

	private boolean exclude;
	private String tag;
	private FilterType filterType;
	// Range variables
	private Double min;
	private Double max;
	private String expression;
	private boolean hasWildCards;
	private boolean isRegExpr;
	private boolean caseSensitive;

	AlignmentFilter() {
		exclude = true;
		tag = "NM";
		filterType = FilterType.missing;
		expression = "*T*";
		hasWildCards = false;
		isRegExpr = false;
		caseSensitive = false;
		min=1.0;
		max=3.0;
	}

	public String getTag() {
		return tag;
	}

	public void setTag(SAMRecord.SAMTagAndValue tag) {
		this.tag = tag.tag;
	}

	public boolean filter(Integer attr) {
		return filter(new Double(attr));
	}

	public boolean filter(int attr) {
		return filter(new Double(attr));
	}

	public boolean filter(Double attr) {
		if (filterType == FilterType.range) {
			if (this.min <= attr && this.max >= attr)
				return true;
			else
				return false;
		} else {
			filter(attr.toString());
		}
		return false;
	}

	public boolean filter(String attr) {
		String regExprToUse = null;
		Pattern regExpr;

		// compile regular expression
		if (hasWildCards) {
			regExprToUse = WildcardMatcher.wildcardToRegex(expression);
		} else if (isRegExpr) {
			regExprToUse = expression;
		} else {
			regExprToUse = null;
		}
		if (regExprToUse != null) {
			// allow - and match - LF and international chars in the data
			int flags = Pattern.DOTALL | Pattern.MULTILINE;
			if (!caseSensitive) {
				flags |= Pattern.CASE_INSENSITIVE | Pattern.UNICODE_CASE;
			}
			regExpr = Pattern.compile(regExprToUse, flags);
		} else {
			regExpr = null;
		}
		boolean match = false;
		if (regExpr != null) {
			Matcher matcher = regExpr.matcher(attr);
			match = matcher.matches();
		}else{
			// otherwise use the string pattern directly
			if (caseSensitive) {
                match = expression.equals(attr);
            } else {
                match = expression.equalsIgnoreCase(attr);
            }
		}

		return match;
	}

	public boolean filter(Object attribute) {
		// TODO Auto-generated method stub
		if (attribute instanceof Integer){
			return filter((Integer) attribute);
		}
		return false;
	}

	//public boolean filter(Object attr) {
	//	return false;
	//}

}
