package org.broad.igv.gwas;

import org.broad.igv.feature.AbstractFeature;
import org.broad.igv.tdf.BufferedByteWriter;
import org.broad.igv.track.WindowFunction;

import java.io.IOException;

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
    public String getValueString(double position, int mouseX, WindowFunction windowFunction) {
        StringBuilder sb = new StringBuilder();
        sb.append(snp);
        sb.append("<br>" + geneId);
        sb.append("<br>" + geneName);
        //sb.append("<br>tStat = " + tStat);
        sb.append("<br>pValue = " + pValue);
        //sb.append("<br>qValue = " + qValue);
        return sb.toString();

    }
}
