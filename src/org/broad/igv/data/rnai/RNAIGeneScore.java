/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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
