package org.broad.igv.gwas;

import org.broad.igv.Globals;
import org.broad.igv.feature.BasicFeature;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;

/**
 * CODEC for GTex project eQTL files
 *
 * SNP	SNP_Chr	SNP_Pos	Gen_ID	Gene_Name	Gene_Pos	T_Stat	P_Val	Q_Val
 * rs1569471	1	169564130	ENSG00000000460.11	C1orf112	169631245	-4.187361378	7.79E-05	0.04794564
 *
 */
public class EQTLCodec  extends AsciiFeatureCodec<EQTLFeature> {


    protected EQTLCodec(Class myClass) {
        super(myClass);
    }

    @Override
    public Feature decodeLoc(String line) {
        String [] tokens = Globals.tabPattern.split(line);
        if(tokens[0].equals("SNP")) return null;
        String chr = tokens[1];
        int position = Integer.parseInt(tokens[2]) - 1;
        return new BasicFeature(chr, position, position+1);
    }

    @Override
    public EQTLFeature decode(String s) {

        String [] tokens = Globals.tabPattern.split(s);
        if(tokens[0].equals("SNP")) return null;

        String snp = tokens[0];
        String chr = tokens[1];
        int position = Integer.parseInt(tokens[2]) - 1;
        String geneId = tokens[3];
        String geneName = tokens[4];
        int genePosition = Integer.parseInt(tokens[5]);
        double tStat = Double.parseDouble(tokens[6]);
        double pValue = Double.parseDouble(tokens[7]);
        double qValue = Double.parseDouble(tokens[8]);

        return new EQTLFeature(snp, chr, position, geneId, geneName, genePosition, tStat, pValue, qValue);
    }
}
