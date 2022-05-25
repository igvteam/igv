package org.broad.igv.bbfile.codecs;

import org.broad.igv.feature.*;
import org.broad.igv.ui.color.ColorUtilities;

import java.awt.*;
import java.util.List;

import static org.broad.igv.feature.tribble.IGVBEDCodec.createExons;
import static org.broad.igv.feature.tribble.IGVBEDCodec.isCoding;

public class BBFeature implements IGVFeature {

    private String chr;
    private int start;
    private int end;


    public BBFeature(String chr, int start, int end){

        this.chr =  chr;
        this.start =  start;
        this.end = end;
        //restOfFields = restOfFieldsString == null ? null : restOfFieldsString.split("\t");
    }

    @Override
    public List<String> getAttributeKeys() {
        return IGVFeature.super.getAttributeKeys();
    }

    @Override
    public String getAttribute(String key) {
        return IGVFeature.super.getAttribute(key);
    }

    @Override
    public List<Exon> getExons() {
        return IGVFeature.super.getExons();
    }

    @Override
    public Color getColor() {
        return null;
    }

    @Override
    public float getScore() {
        return 0;
    }

    @Override
    public String getName() {
        return null;
    }

    @Override
    public String getContig() {
        return null;
    }

    @Override
    public int getStart() {
        return 0;
    }

    @Override
    public int getEnd() {
        return 0;
    }

/*    public void decodeRestOfFiels() {

        // The first 3 columns are non optional for BED.  We will relax this
        // and only require 2.
        int tokenCount = tokens.length;

        if (tokenCount < 2) {
            return null;
        }

        String c = tokens[0];
        String chr = genome == null ? c : genome.getCanonicalChrName(c);

        //BED format, and IGV, use starting element as 0.
        int start = Integer.parseInt(tokens[1]);

        int end = start + 1;
        if (tokenCount > 2) {
            end = Integer.parseInt(tokens[2]);
        }

        BasicFeature feature = new BasicFeature(chr, start, end);

        // The rest of the columns are optional.  Stop parsing upon encountering
        // a non-expected value

        // Name
        if (tokenCount > 3) {
            String name = tokens[3].replaceAll("\"", "");
            if (name.equals(".")) name = "";   // Convention
            feature.setName(name);
            feature.setIdentifier(name);
        }


        // Bed files are not always to-spec after the name field.  Stop parsing when we find an unexpected column.
        // Score

        if (tokenCount > 4) {
            try {
                float score = tokens[4].equals(".") ? 1000 : Float.parseFloat(tokens[4]);
                feature.setScore(score);
            } catch (NumberFormatException numberFormatException) {
                // Interpret as missing score value, set to 1000 => max value for score according to BED spec.
                feature.setScore(1000);
            }
        }

        // Strand
        if (tokenCount > 5) {
            String strandString = tokens[5].trim();
            char strand = (strandString.length() == 0)
                    ? ' ' : strandString.charAt(0);

            if (strand == '-') {
                feature.setStrand(Strand.NEGATIVE);
            } else if (strand == '+') {
                feature.setStrand(Strand.POSITIVE);
            } else {
                feature.setStrand(Strand.NONE);
            }
        }

        // Thick ends
        if (tokenCount > 7) {
            try {
                int thickStart = Integer.parseInt(tokens[6]);
                int thickEnd = Integer.parseInt(tokens[7]);
                feature.setThickStart(Math.max(start, thickStart));
                feature.setThickEnd(Math.min(end, thickEnd));
            } catch (NumberFormatException e) {
                return feature;
            }
        }


        // Color
        if (tokenCount > 8) {
            String colorString = tokens[8];
            if (colorString.trim().length() > 0 && !colorString.equals(".")) {
                feature.setColor(ColorUtilities.stringToColor(colorString));
            }
        }

        // Exons
        if (tokenCount > 11) {
            createExons(start, tokens, feature, chr, feature.getStrand());
        }

        if(isCoding(feature)) {
            FeatureUtils.computeReadingFrames(feature);
        }

        return feature;
    }*/

}
