package org.igv.circview.util;

import java.awt.Color;

/**
 * Parsing and manipulation of CSS-style color strings used by the JS source
 * ({@code rgb(...)}, {@code rgba(...)}, and {@code #hex}).
 *
 * <p>The alpha helpers mirror setAlpha()/getAlpha() in circularView.js, where
 * alpha is a fraction in [0, 1].
 */
public final class ColorUtils {

    private ColorUtils() {
    }

    /**
     * Parse a CSS color string into an AWT Color.
     * Supports {@code rgb(r,g,b)}, {@code rgba(r,g,b,a)} (a in [0,1]),
     * {@code #rgb}, and {@code #rrggbb}, plus a handful of common named colors.
     *
     * @throws IllegalArgumentException if the string cannot be parsed
     */
    public static Color parseColor(String s) {
        if (s == null) {
            throw new IllegalArgumentException("color string is null");
        }
        String c = s.trim();
        String lower = c.toLowerCase();

        if (lower.startsWith("rgba(") || lower.startsWith("rgb(")) {
            int open = c.indexOf('(');
            int close = c.indexOf(')');
            String body = close > open ? c.substring(open + 1, close) : c.substring(open + 1);
            String[] parts = body.split(",");
            int r = clamp255((int) Math.round(parseNum(parts[0])));
            int g = clamp255((int) Math.round(parseNum(parts[1])));
            int b = clamp255((int) Math.round(parseNum(parts[2])));
            int a = 255;
            if (parts.length >= 4) {
                a = clamp255((int) Math.round(parseNum(parts[3]) * 255.0));
            }
            return new Color(r, g, b, a);
        }

        if (c.startsWith("#")) {
            String hex = c.substring(1);
            if (hex.length() == 3) {
                int r = Integer.parseInt(hex.substring(0, 1), 16) * 17;
                int g = Integer.parseInt(hex.substring(1, 2), 16) * 17;
                int b = Integer.parseInt(hex.substring(2, 3), 16) * 17;
                return new Color(r, g, b);
            } else if (hex.length() == 6) {
                int r = Integer.parseInt(hex.substring(0, 2), 16);
                int g = Integer.parseInt(hex.substring(2, 4), 16);
                int b = Integer.parseInt(hex.substring(4, 6), 16);
                return new Color(r, g, b);
            }
        }

        switch (lower) {
            case "black":
                return Color.BLACK;
            case "white":
                return Color.WHITE;
            case "red":
                return Color.RED;
            case "green":
                return Color.GREEN;
            case "blue":
                return Color.BLUE;
            case "gray":
            case "grey":
                return Color.GRAY;
            default:
                throw new IllegalArgumentException("Unrecognized color: " + s);
        }
    }

    /** Return a copy of {@code color} with the given alpha fraction in [0, 1]. */
    public static Color setAlpha(Color color, float alpha) {
        int a = clamp255(Math.round(alpha * 255f));
        return new Color(color.getRed(), color.getGreen(), color.getBlue(), a);
    }

    /** Alpha of {@code color} as a fraction in [0, 1]. */
    public static float getAlpha(Color color) {
        return color.getAlpha() / 255f;
    }

    private static double parseNum(String s) {
        return Double.parseDouble(s.trim());
    }

    private static int clamp255(int v) {
        return Math.max(0, Math.min(255, v));
    }
}
