package org.broad.igv.gwas;

import org.broad.tribble.Feature;

/**
 * Represents an eQTL value
 * <p/>
 * SNP	SNP_Chr	SNP_Pos	Gen_ID	Gene_Name	Gene_Pos	T_Stat	P_Val	Q_Val
 * rs1569471	1	169564130	ENSG00000000460.11	C1orf112	169631245	-4.187361378	7.79E-05	0.04794564
 */
public class EQTLFeature implements Feature {

    String snp;
    String chr;
    int position;
    String geneId;
    String geneName;
    int genePosition;
    double tStat;
    double pValue;
    double qValue;

    public EQTLFeature(String snp, String chr, int position, String geneId, String geneName,
                       int genePosition, double tStat, double pValue, double qValue) {
        this.snp = snp;
        this.chr = chr;
        this.position = position;
        this.geneId = geneId;
        this.geneName = geneName;
        this.genePosition = genePosition;
        this.tStat = tStat;
        this.pValue = pValue;
        this.qValue = qValue;
    }

    @Override
    public String getChr() {
        return chr;
    }

    @Override
    public int getStart() {
        return position;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public int getEnd() {
        return position + 1;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
