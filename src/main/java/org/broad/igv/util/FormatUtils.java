package org.broad.igv.util;

import java.util.List;
import java.util.Map;

public class FormatUtils {

    private static final int MAX_CHARS_PER_LINE = 50;

    public static void printHtml(Map<String, String> map, StringBuffer buffer, int max) {

        if (map == null || map.isEmpty()) return;

        int count = 0;
        buffer.append("<br>");
        for (Map.Entry<String, String> entry : map.entrySet()) {
            String value = entry.getValue();

            buffer.append("<b>" + entry.getKey() + "</b>");
            buffer.append(":&nbsp;");
            String ts;
            if (value.startsWith("https://")) {
                ts = "<a href='" + value + "'>" + value + "</a>";
            } else {
                ts = value.startsWith("<") ? value : lineWrapString(value, MAX_CHARS_PER_LINE);
            }

            buffer.append(ts);
            buffer.append("<br/>");
            count++;

            if (++count > max) {
                buffer.append("...");
                break;
            }

        }
    }

    private static String lineWrapString(String input, int maxCharsPerLine) {
        int lines = input.length() / maxCharsPerLine + 1;
        if (lines == 1) return input;

        String result = input.substring(0, maxCharsPerLine);
        for (int lineNum = 1; lineNum < lines; lineNum++) {
            int start = lineNum * maxCharsPerLine;
            int end = Math.min(start + maxCharsPerLine, input.length());
            result += "<br/>" + input.substring(start, end);
        }
        return result;
    }


}
