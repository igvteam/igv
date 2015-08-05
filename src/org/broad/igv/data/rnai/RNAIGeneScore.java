/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */


package org.broad.igv.data.rnai;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.NamedFeature;
import org.broad.igv.track.WindowFunction;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Class description
 *
 * @author Enter your name here...
 * @version Enter version here..., 08/10/06
 */
public class RNAIGeneScore implements LocusScore {

    String batchCondition;
    NamedFeature gene;
    float geneScore;
    int geneConfidence;
    boolean hasGeneConfidence;
    int numberOfHairpins;
    int start;
    int end;

    /**
     * Constructs ...
     *
     * @param batchId
     * @param gene
     * @param geneScore
     * @param numberOfHairpins
     */
    public RNAIGeneScore(String batchId, NamedFeature gene, float geneScore,
                         int numberOfHairpins) {
        this.batchCondition = batchId;
        this.gene = gene;
        this.geneScore = geneScore;
        geneConfidence = 3;
        hasGeneConfidence = false;
        this.numberOfHairpins = numberOfHairpins;
        this.start = gene.getStart();
        this.end = gene.getEnd();
    }

    /**
     * Constructs ...
     *
     * @param batchId
     * @param gene
     * @param geneScore
     * @param confidence
     * @param numberOfHairpins
     */
    public RNAIGeneScore(String batchId, NamedFeature gene, float geneScore, int confidence,
                         int numberOfHairpins) {
        this.batchCondition = batchId;
        this.gene = gene;
        this.geneScore = geneScore;
        this.geneConfidence = confidence;
        hasGeneConfidence = true;
        this.numberOfHairpins = numberOfHairpins;
        this.start = gene.getStart();
        this.end = gene.getEnd();
    }

    /**
     * Copy constructor
     *
     * @param score
     */
    public RNAIGeneScore(RNAIGeneScore score) {
        this.batchCondition = score.batchCondition;
        this.gene = score.gene;
        this.geneScore = score.geneScore;
        this.geneConfidence = score.geneConfidence;
        this.numberOfHairpins = score.numberOfHairpins;
        this.start = score.getStart();
        this.end = score.getEnd();
    }

    /**
     * Method description
     *
     * @return
     */
    public RNAIGeneScore copy() {
        return new RNAIGeneScore(this);
    }


    /**
     * Checks if there is a confidence value for the gene score
     *
     * @return
     */
    public boolean hasGeneConfidence() {
        return hasGeneConfidence;
    }

    /**
     * Method description
     *
     * @return
     */
    public NamedFeature getGene() {
        return gene;
    }

    public String getChr() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getContig() {
        return null;
    }

    /**
     * Method description
     *
     * @return
     */
    public int getStart() {
        return start;
    }

    /**
     * Method description
     *
     * @return
     */
    public int getEnd() {
        return end;
    }

    /**
     * Method description
     *
     * @return
     */
    public float getScore() {
        return geneScore;
    }

    public float getConfidence() {
        return geneConfidence;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    /**
     * Return the list of hairpin scores associated with this gene.  Scores are sorted alphabetically
     * by haipin name
     *
     * @return
     */
    List<RNAIHairpinValue> hairpinScores = null;

    /**
     * Method description
     *
     * @return
     */
    public Collection<RNAIHairpinValue> getHairpinValues() {
        if (hairpinScores == null) {
            Collection<RNAIHairpinValue> tmp =
                    RNAIHairpinCache.getInstance().getHairpinScores(batchCondition, gene.getName().toUpperCase());
            if (tmp != null) {
                hairpinScores = new ArrayList(tmp);
            }
        }

        return hairpinScores;
    }

    /**
     * Method description
     *
     * @param ignored
     * @return
     */
    public String getValueString(double position, WindowFunction ignored) {

        StringBuffer buf = new StringBuffer(100);
        buf.append("<b>Gene: " + gene.getName() + "</b><br>");
        buf.append("<b>Intensity: " + geneScore + "<br></b>");

        if (hasGeneConfidence)
            buf.append("Confidence: " + geneConfidence + "<br>");

        Collection<RNAIHairpinValue> hpins = getHairpinValues();
        if ((hpins == null) || hpins.isEmpty()) {
            buf.append("# Hairpins: " + numberOfHairpins + "<br>");
        } else {
            buf.append("Hairpin scores:<br>");
            for (RNAIHairpinValue hpScore : hpins) {
                if (hpScore.hasScoreSTD()) {
                    buf.append(hpScore.getName() + ": " + hpScore.getScoreMean() + "  ("
                            + hpScore.getScoreSTD() + ")<br>");
                } else {
                    buf.append(hpScore.getName() + ": " + hpScore.getScoreMean() + "<br>");
                }
            }
        }

        return buf.toString();

    }
}
