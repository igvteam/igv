package org.igv.variant;

import htsjdk.variant.vcf.VCFHeader;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.renderer.ColorStopScale;
import org.igv.variant.vcf.VCFVariant;

import java.awt.*;
import java.util.List;

/**
 * Coloring variants by allele frequency so that rare variants stand out: the top-level "Allele Frequency" (AF or
 * GMAF) and "Allele Fraction" (AC / AN) modes, and well-known allele frequency INFO fields such as gnomAD_AF.
 * Frequencies span orders of magnitude, so colors come from a {@link ColorStopScale}, blended on a log scale between
 * stops.
 */
public class AlleleFrequencyColors {

    private static Logger log = LogManager.getLogger(AlleleFrequencyColors.class);

    /**
     * Values of the {@link Constants#VARIANT_ALLELE_FREQUENCY_DISPLAY} preference.
     */
    public static final String DISPLAY_AUTOMATIC = "AUTOMATIC";
    public static final String DISPLAY_COLOR_SCALE = "COLOR_SCALE";
    public static final String DISPLAY_BAR = "BAR";

    /**
     * Well-known allele frequency INFO fields, colored by rarity by default when chosen under color by INFO field.
     * AF is not among them -- it has the top-level "Allele Frequency" item.
     */
    public static final List<String> FREQUENCY_INFO_FIELDS = List.of(
            "GMAF", "gnomAD_AF", "gnomADe_AF", "gnomADg_AF", "gnomAD_exomes_AF", "gnomAD_genomes_AF",
            "AF_joint", "AF_grpmax", "AF_popmax", "MAX_AF", "ExAC_AF");

    /**
     * Default stops, following gnomAD usage.  Common variants (5% and above) are a neutral blue-gray; low frequency
     * (1-5%) shifts toward yellow; rare (0.1-1%) to orange; 0.01% is red; 0.0001% magenta -- about the frequency of
     * a single allele among gnomAD v4's roughly 1.6 million.
     */
    public static final ColorStopScale DEFAULT_SCALE = new ColorStopScale(List.of(
            new ColorStopScale.Stop(1, new Color(130, 145, 175)),
            new ColorStopScale.Stop(0.05, new Color(130, 145, 175)),
            new ColorStopScale.Stop(0.01, new Color(240, 200, 40)),
            new ColorStopScale.Stop(0.001, new Color(245, 130, 30)),
            new ColorStopScale.Stop(0.0001, new Color(220, 30, 30)),
            new ColorStopScale.Stop(0.000001, new Color(200, 0, 200))));

    private static String cachedScaleString;
    private static ColorStopScale cachedScale;

    /**
     * @return the scale set by the user, or the default
     */
    public static synchronized ColorStopScale getScale() {
        String value = PreferencesManager.getPreferences().get(Constants.VARIANT_ALLELE_FREQUENCY_COLORS, null);
        if (value == null || value.isBlank()) {
            return DEFAULT_SCALE;
        }
        if (!value.equals(cachedScaleString)) {
            try {
                cachedScale = new ColorStopScale(value);
            } catch (RuntimeException e) {
                log.error("Invalid allele frequency colors preference: " + value, e);
                cachedScale = DEFAULT_SCALE;
            }
            cachedScaleString = value;
        }
        return cachedScale;
    }

    /**
     * Save the scale as a preference, applying to all variant tracks.  The default scale is saved by removing the
     * preference, so that a user who never changed the colors gets any future change to the defaults.
     */
    public static synchronized void setScale(ColorStopScale scale) {
        if (scale == null || scale.asString().equals(DEFAULT_SCALE.asString())) {
            PreferencesManager.getPreferences().remove(Constants.VARIANT_ALLELE_FREQUENCY_COLORS);
        } else {
            PreferencesManager.getPreferences().put(Constants.VARIANT_ALLELE_FREQUENCY_COLORS, scale.asString());
        }
    }

    /**
     * @return true if the INFO field is a well-known allele frequency field, colored by rarity by default
     */
    public static boolean isFrequencyField(String infoKey) {
        return infoKey != null && FREQUENCY_INFO_FIELDS.stream().anyMatch(infoKey::equalsIgnoreCase);
    }

    /**
     * The frequency to color a variant by in the top-level modes: for ALLELE_FREQUENCY the AF value, or GMAF if there
     * is no AF -- the fields the bar display reads -- and for ALLELE_FRACTION AC / AN.  At a multi-allelic site the
     * rarest alternate allele is used.
     *
     * @return the frequency, or -1 if the variant has no value
     */
    public static double getFrequency(Variant variant, VariantTrack.ColorMode mode) {

        if (mode == VariantTrack.ColorMode.ALLELE_FRACTION) {
            int[] counts = variant.getAlleleCounts();
            int total = variant.getTotalAlleleCount();
            if (counts == null || counts.length == 0 || total <= 0) {
                return -1;
            }
            int min = counts[0];
            for (int count : counts) {
                min = Math.min(min, count);
            }
            return (double) min / total;
        }

        for (String key : VCFVariant.ALLELE_FREQUENCY_KEYS) {
            double frequency = minimumValue(variant.getAttributeAsString(key));
            if (frequency >= 0) {
                return frequency;
            }
        }
        return -1;
    }

    /**
     * @return true if the header declares the field(s) the mode takes values from.  If it doesn't, a missing value
     * means the file has no frequencies at all, not that the variant is unobserved, so it shouldn't be colored as rare.
     */
    public static boolean hasFrequencyFields(VCFHeader header, VariantTrack.ColorMode mode) {
        if (mode == VariantTrack.ColorMode.ALLELE_FRACTION) {
            return header.getInfoHeaderLine("AC") != null && header.getInfoHeaderLine("AN") != null;
        }
        for (String key : VCFVariant.ALLELE_FREQUENCY_KEYS) {
            if (header.getInfoHeaderLine(key) != null) {
                return true;
            }
        }
        return false;
    }

    /**
     * The smallest non-negative number in a comma separated list, ignoring missing (".") entries.  Brackets are
     * stripped, as htsjdk may present a list attribute as "[0.1, 0.2]".
     *
     * @return the minimum, or -1 if there is none
     */
    private static double minimumValue(String value) {
        if (value == null) {
            return -1;
        }
        double min = -1;
        for (String token : value.replaceAll("[\\[\\]\\(\\)]", "").split(",")) {
            try {
                double number = Double.parseDouble(token.trim());
                if (number >= 0 && (min < 0 || number < min)) {
                    min = number;
                }
            } catch (NumberFormatException e) {
                // "." -- no value for this allele
            }
        }
        return min;
    }
}
