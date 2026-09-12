package org.igv.variant;

import org.igv.ui.color.ColorUtilities;
import org.igv.ui.color.PaletteColorTable;

import java.awt.*;

/**
 * Color tables used when coloring variants by a VCF INFO attribute.
 * <p>
 * A few widely used attributes (SVTYPE, VT, CLNSIG) have predefined tables so that, for example, a deletion is
 * always red and a pathogenic ClinVar variant is always red.  Any other attribute gets colors from a categorical
 * palette, assigned as values are encountered.  Values not covered by a predefined table also fall back to the
 * palette, so a table can be incomplete without breaking anything.
 */
class VariantColorTables {

    /**
     * Palette used for attributes with no predefined table.  Matches the default used by igv.js.
     */
    private static final String DEFAULT_PALETTE = "Set 1";

    private VariantColorTables() {
    }

    /**
     * Return a new (mutable) color table for the given INFO attribute.  Each track gets its own instance, as
     * tables accumulate color assignments as values are encountered.
     */
    static PaletteColorTable getColorTable(String infoKey) {
        switch (infoKey == null ? "" : infoKey.toUpperCase()) {
            case "SVTYPE":
                return svTypeTable();
            case "VT":
                return variantTypeTable();
            case "CLNSIG":
                return clinicalSignificanceTable();
            default:
                return newPaletteTable();
        }
    }

    private static PaletteColorTable newPaletteTable() {
        return new PaletteColorTable(ColorUtilities.getPalette(DEFAULT_PALETTE));
    }

    /**
     * Structural variant types.  Colors match igv.js.
     */
    private static PaletteColorTable svTypeTable() {
        PaletteColorTable table = newPaletteTable();
        table.put("DEL", new Color(255, 33, 1));
        table.put("INS", new Color(0, 24, 136));
        table.put("DUP", new Color(2, 132, 1));
        table.put("INV", new Color(0, 134, 136));
        table.put("CNV", new Color(137, 49, 255));
        table.put("BND", new Color(137, 17, 0));
        return table;
    }

    /**
     * Variant type as annotated by the 1000 Genomes project.
     */
    private static PaletteColorTable variantTypeTable() {
        PaletteColorTable table = newPaletteTable();
        table.put("SNP", new Color(55, 126, 184));
        table.put("INDEL", new Color(228, 26, 28));
        table.put("SV", new Color(77, 175, 74));
        table.put("MNP", new Color(152, 78, 163));
        return table;
    }

    /**
     * ClinVar clinical significance, ordered benign (green) to pathogenic (red), with uncertain and conflicting
     * classifications in neutral colors.  Both the current ("classifications") and legacy ("interpretations")
     * spellings of the conflicting value are included.
     */
    private static PaletteColorTable clinicalSignificanceTable() {

        final Color pathogenic = new Color(202, 0, 32);
        final Color likelyPathogenic = new Color(244, 109, 67);
        final Color uncertain = new Color(150, 150, 150);
        final Color conflicting = new Color(230, 171, 2);
        final Color likelyBenign = new Color(146, 197, 222);
        final Color benign = new Color(5, 113, 176);

        PaletteColorTable table = newPaletteTable();
        table.put("Pathogenic", pathogenic);
        table.put("Pathogenic/Likely_pathogenic", pathogenic);
        table.put("Likely_pathogenic", likelyPathogenic);
        table.put("Uncertain_significance", uncertain);
        table.put("Conflicting_classifications_of_pathogenicity", conflicting);
        table.put("Conflicting_interpretations_of_pathogenicity", conflicting);
        table.put("Likely_benign", likelyBenign);
        table.put("Benign/Likely_benign", benign);
        table.put("Benign", benign);
        return table;
    }
}
