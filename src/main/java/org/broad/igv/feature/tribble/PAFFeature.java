package org.broad.igv.feature.tribble;

import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.Strand;
import org.broad.igv.track.WindowFunction;

import java.awt.*;

/**
 * Created by jrobinso on 10/31/17.
 */
public class PAFFeature implements IGVFeature {

    static Color defaultColor = new Color(0, 0, 150);

    private static int chrColumn = 5;
    private static int startColumn = 7;
    private static int endColumn = 8;
    private static int scoreColumn = 11;


    String chr;
    int start;
    int end;
    Strand strand;
    float score;
    String name;
    String description;


    public PAFFeature(String chr, int start, int end, Strand strand, float score, String name, String description) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.score = score;
        this.name = name;
        this.description = description;
    }


    @Override
    public String getType() {
        return "Aignment";
    }

    @Override
    public String getIdentifier() {
        return name;
    }

    @Override
    public String getDescription() {
        return null;
    }

    @Override
    public Strand getStrand() {
        return strand;
    }

    @Override
    public Color getColor() {
        return defaultColor;
    }

    @Override
    public void setStart(int start) {
        this.start = start;
    }

    @Override
    public void setEnd(int end) {
        this.end = end;
    }

    @Override
    public float getScore() {
        return score;
    }

    @Override
    public String getValueString(double position, int mouseX, WindowFunction windowFunction) {
        return description;
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public String getContig() {
        return chr;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getEnd() {
        return end;
    }
}
