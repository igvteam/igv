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

package org.broad.igv.gwas;

import org.broad.igv.feature.AbstractFeature;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.Strand;
import org.broad.igv.tdf.BufferedByteWriter;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.collections.MultiMap;
import htsjdk.tribble.Feature;

import java.awt.*;
import java.io.IOException;
import java.util.List;
import java.util.Map;

/**
 * Represents an eQTL value
 * <p/>
 * SNP	SNP_Chr	SNP_Pos	Gen_ID	Gene_Name	Gene_Pos	T_Stat	P_Val	Q_Val
 * rs1569471	1	169564130	ENSG00000000460.11	C1orf112	169631245	-4.187361378	7.79E-05	0.04794564
 *
 * SNP	SNP_Chr	SNP_Pos	Gen_ID	Gene_Name	Gene_Pos	T_Stat	Beta	P_Val	min(p)	EmpP	nom_thresh
 chr2:202672143:I	2	202672143	ENSG00000003393.10	ALS2	202645912	-4.90641599377695	-0.410144871510459	4.16077158322729e-06	4.1199296E-8	9 *

 */
public class EQTLFeature extends AbstractFeature {

    private String snp;
    String chr;
    int position;
    private String geneId;
    private String geneName;
    //private float tStat;
    //private float beta;
    private float pValue;
    //private float qValue;

    // function Eqtl(snp, chr, position, geneId, geneName, pValue)
    public EQTLFeature(String snp, String chr, int position, String geneId, String geneName, float pValue) {
        this.snp = snp;
        this.chr = chr;
        this.position = position;
        this.geneId = geneId;
        this.geneName = geneName;
        //this.tStat = tStat;
        this.pValue = pValue;
        //this.qValue = qValue;
    }


    public byte [] encodeBinary() throws IOException {
        BufferedByteWriter writer = new BufferedByteWriter();
        writer.putNullTerminatedString(snp);
        writer.putNullTerminatedString(chr);
        writer.putInt(position);
        writer.putNullTerminatedString(geneId);
        writer.putNullTerminatedString(geneName);
        //writer.putFloat(tStat);
        writer.putFloat(pValue);
        //writer.putFloat(qValue);
        return writer.getBytes();
    }

    @Override
    public String getChr() {
        return chr;
    }

    @Override
    public int getStart() {
        return position;
    }

    @Override
    public int getEnd() {
        return position + 1;
    }


    public String getSnp() {
        return snp;
    }

    public String getGeneId() {
        return geneId;
    }

    public String getGeneName() {
        return geneName;
    }

    @Override
    public String getValueString(double position, WindowFunction windowFunction) {
        StringBuilder sb = new StringBuilder();
        sb.append(snp);
        sb.append("<br>" + geneId);
        sb.append("<br>" + geneName);
        //sb.append("<br>tStat = " + tStat);
        sb.append("<br>pValue = " + pValue);
        //sb.append("<br>qValue = " + qValue);
        return sb.toString();

    }

    @Override
    public String getURL() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
