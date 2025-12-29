package org.igv.ucsc.bb.codecs;

import org.igv.Globals;
import org.igv.ucsc.bb.BBUtils;

import org.igv.feature.BasicFeature;
import org.igv.feature.Exon;
import org.igv.feature.Strand;
import org.igv.ucsc.bb.BedData;

import java.awt.*;
import java.util.Map;

public class BBRmskCodec implements BBCodec {

    BBBedCodec bedCodec;

    public BBRmskCodec(int standardFieldCount, BBUtils.ASTable autosql) {
        bedCodec = new BBBedCodec(standardFieldCount, autosql);
    }

    @Override
    public BasicFeature decode(BedData data) {

        BasicFeature feature = bedCodec.decode(data);

        String id = feature.getAttribute("id");
        if (id != null) {
            feature.setIdentifier(id);
            feature.removeAttribute("id");
        }

        String name = feature.getName();
        String repClass = null;
        String[] parts = name.split("#");
        if (parts.length > 1) {
            name = parts[0];
            repClass = parts[1];
        }
        feature.setName(name);

        if (repClass != null) {
            feature.setAttribute("Repeat Class", repClass);
            int idx1 = repClass.indexOf("/");
            String c = idx1 > 0 ? repClass.substring(0, idx1) : repClass;
            if (c.endsWith("?")) c = c.substring(c.length() - 1);
            Color color = classColors.containsKey(c) ? classColors.get(c) : defaultColor;
            feature.setColor(color);
        }


        int blockCount = Integer.parseInt(feature.getAttribute("blockCount"));
        createExons(feature, blockCount,
                Globals.commaPattern.split(feature.getAttribute("blockSizes")),
                Globals.commaPattern.split(feature.getAttribute("blockStarts")));

        feature.removeAttribute("blockCount");
        feature.removeAttribute("blockSizes");
        feature.removeAttribute("blockStarts");


        // Compatibility with UCSC,  meaning of "start" is ambiguous, but UCSC visualization draws from "thickStart"
        feature.setStart(feature.getThickStart());
        feature.setEnd(feature.getThickEnd());


        return feature;
    }

    /**
     * Compute blocks, basically a copy of the equivalent method in IGVBedCodec.  The bigRmsk format includes
     * the block information in non-standard fields, presumably because "itemRGB" is skipped.
     */
    private static void createExons(BasicFeature feature, int blockCount, String[] blockSizes, String[] blockStarts) throws NumberFormatException {

        int start = feature.getStart();
        String chr = feature.getChr();
        Strand strand = feature.getStrand();

        if (blockStarts.length == blockSizes.length) {
            for (int i = 0; i < blockStarts.length; i++) {
                int bstart = Integer.parseInt(blockStarts[i]);
                int bsize = Integer.parseInt(blockSizes[i]);
                if (bstart >= 0 && bsize > 0) {
                    int exonStart = start + bstart;
                    int exonEnd = exonStart + bsize;
                    if(exonStart < feature.getThickStart()) {
                        feature.setThickStart(exonStart);
                    }
                    if(exonEnd > feature.getThickEnd()) {
                        feature.setThickEnd(exonEnd);
                    }
                    Exon exon = new Exon(chr, exonStart, exonEnd, strand);
                    feature.addExon(exon);
                }
            }
        }
    }



    static Map<String, Color> classColors = Map.of(
            "SINE", new Color(31, 119, 180),
            "LINE", new Color(255, 127, 14),
            "LTR", new Color(44, 160, 44),
            "DNA", new Color(214, 39, 40),
            "Simple_repeat", new Color(148, 103, 189),
            "Low_complexity", new Color(140, 86, 75),
            "Satellite", new Color(227, 119, 194),
            "RNA", new Color(127, 127, 127),
            "Other", new Color(188, 189, 34),
            "Unknown", new Color(23, 190, 207)
    );

    static Color defaultColor = new Color(0, 180, 180);
}
