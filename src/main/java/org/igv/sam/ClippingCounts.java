package org.igv.sam;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class ClippingCounts {
    private final static ClippingCounts NO_CLIPPING = new ClippingCounts(0, 0, 0, 0);
    public static final Pattern LCLIP_PATTERN = Pattern.compile("^(([0-9]+)H)?(([0-9]+)S)?");
    public static final Pattern RCLIP_PATTERN = Pattern.compile("(([0-9]+)S)?(([0-9]+)H)?$");
    private final int leftHard;
    private final int leftSoft;
    private final int rightHard;
    private final int rightSoft;

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (!(o instanceof ClippingCounts)) return false;

        final ClippingCounts that = (ClippingCounts) o;

        if (leftHard != that.leftHard) return false;
        if (leftSoft != that.leftSoft) return false;
        if (rightHard != that.rightHard) return false;
        return rightSoft == that.rightSoft;
    }

    @Override
    public int hashCode() {
        int result = leftHard;
        result = 31 * result + leftSoft;
        result = 31 * result + rightHard;
        result = 31 * result + rightSoft;
        return result;
    }

    public ClippingCounts(final int leftHard, final int leftSoft, final int rightSoft, final int rightHard) {
        this.leftHard = leftHard;
        this.leftSoft = leftSoft;
        this.rightHard = rightHard;
        this.rightSoft = rightSoft;
    }

    public static ClippingCounts fromCigar(Cigar cigar) {
        if (cigar == null || cigar.isEmpty() || !cigar.isClipped()) {
            return NO_CLIPPING;
        } else {
            final int numberOfElements = cigar.numCigarElements();
            int leftHard = 0, leftSoft = 0, rightSoft = 0, rightHard = 0;
            int i = 0;
            for (; i < 2 && i <= numberOfElements; i++) {
                final CigarElement element = cigar.getCigarElement(i);
                final CigarOperator operator = element.getOperator();
                if (operator == CigarOperator.H) {
                    leftHard = element.getLength();
                } else if (operator == CigarOperator.S) {
                    leftSoft = element.getLength();
                    break;
                } else {
                    break;
                }
            }
            for (int j = numberOfElements - 1; j > (numberOfElements - 3) && (j > i) && (j > 0); j--) {
                final CigarElement element = cigar.getCigarElement(j);
                final CigarOperator operator = element.getOperator();
                if (operator == CigarOperator.H) {
                    rightHard = element.getLength();
                } else if (operator == CigarOperator.S) {
                    rightSoft = element.getLength();
                    break;
                } else {
                    break;
                }
            }
            return new ClippingCounts(leftHard, leftSoft, rightSoft, rightHard);
        }
    }

    public static ClippingCounts fromCigarString(String cigarString) {
        // Identify the number of hard and soft clipped bases.
        Matcher lclipMatcher = LCLIP_PATTERN.matcher(cigarString);
        Matcher rclipMatcher = RCLIP_PATTERN.matcher(cigarString);
        int lclipHard = 0, lclipSoft = 0, rclipHard = 0, rclipSoft = 0;
        if (lclipMatcher.find()) {
            lclipHard = lclipMatcher.group(2) == null ? 0 : Integer.parseInt(lclipMatcher.group(2), 10);
            lclipSoft = lclipMatcher.group(4) == null ? 0 : Integer.parseInt(lclipMatcher.group(4), 10);
        }
        if (rclipMatcher.find()) {
            rclipHard = rclipMatcher.group(4) == null ? 0 : Integer.parseInt(rclipMatcher.group(4), 10);
            rclipSoft = rclipMatcher.group(2) == null ? 0 : Integer.parseInt(rclipMatcher.group(2), 10);
        }
        return new ClippingCounts(lclipHard, lclipSoft, rclipSoft, rclipHard);
    }

    public boolean isClipped() {
        return leftHard + leftSoft + rightHard + rightSoft > 0;
    }

    public int getLeftHard() {
        return leftHard;
    }

    public int getLeftSoft() {
        return leftSoft;
    }

    public int getRightHard() {
        return rightHard;
    }

    public int getRightSoft() {
        return rightSoft;
    }

    public boolean isLeftClipped() {
        return getLeft() > 0;
    }

    public int getLeft() {
        return leftSoft + leftHard;
    }

    public boolean isRightClipped() {
        return getRight() > 0;
    }

    public int getRight() {
        return rightSoft + rightHard;
    }

    @Override
    public String toString() {
        return "<" + leftHard + "H" + leftSoft + "S" + "|" + rightSoft + "S" + rightHard + "H" + ">";
    }
}
