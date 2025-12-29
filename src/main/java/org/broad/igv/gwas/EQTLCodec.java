package org.broad.igv.gwas;

import org.broad.igv.Globals;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.FeatureFileHeader;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.util.StringUtils;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.readers.LineIterator;

import java.util.HashMap;
import java.util.Map;

/**
 * CODEC for GTex project eQTL files
 * <p/>
 * SNP	SNP_Chr	SNP_Pos	Gen_ID	Gene_Name	Gene_Pos	T_Stat	P_Val	Q_Val
 * rs1569471	1	169564130	ENSG00000000460.11	C1orf112	169631245	-4.187361378	7.79E-05	0.04794564
  SNP	SNP_Chr	SNP_Pos	Gen_ID	Gene_Name	Gene_Pos	T_Stat	Beta	P_Val	min(p)	EmpP	nom_thresh
 chr2:202672143:I	2	202672143	ENSG00000003393.10	ALS2	202645912	-4.90641599377695	-0.410144871510459	4.16077158322729e-06	4.1199296E-8	9 *

 */
public class EQTLCodec extends AsciiFeatureCodec<EQTLFeature> {

    String[] columnNames;
    Genome genome;
    private FeatureFileHeader header;

    public EQTLCodec(Genome genome) {
        super(EQTLFeature.class);
        this.genome = genome;
    }


    //@Override
    public Feature decodeLoc(String line) {
        String[] tokens = Globals.tabPattern.split(line);
        if (tokens[0].equals("SNP")) return null;
        String chr = tokens[1];
        int position = Integer.parseInt(tokens[2]) - 1;
        return new BasicFeature(chr, position, position + 1);
    }

    /*
            var snp = parser.getString();
        var chr = parser.getString();
        var position = parser.getInt();
        var geneId = parser.getString();
        var geneName = parser.getString();
        //var genePosition = -1;
        //var fStat = parser.getFloat();
        var pValue = parser.getFloat();

     */

    @Override
    public EQTLFeature decode(String s) {

        String[] tokens = Globals.tabPattern.split(s);
        if (tokens[0].equals("SNP")) {
            // This is the header
            columnNames = tokens;
            return null;
        }
//  SNP	SNP_Chr	SNP_Pos	Gen_ID	Gene_Name	Gene_Pos	T_Stat	Beta	P_Val	min(p)	EmpP	nom_thresh

        String snp = tokens[0];
        String chr = genome == null ? tokens[1] : genome.getCanonicalChrName(tokens[1]);

        int position = Integer.parseInt(tokens[2]) - 1;
        String geneId = tokens[3];
        String geneName = tokens[4];

        //float tStat = Float.parseFloat(tokens[6]);
        //float beta = Float.parseFloat(tokens[6]);

        double tmp = Double.parseDouble(tokens[8]);
        float pValue = tmp < Float.MIN_VALUE ? Float.MIN_VALUE : (float) tmp;


        //float qValue = tokens.length > 8 ? Float.parseFloat(tokens[8]) : 0f;

        Map<String, String> attributes = null;
        if (columnNames != null) {
            attributes = new HashMap<String, String>();
            for (int i = 5; i < tokens.length; i++) {
                if (columnNames.length < i) {
                    attributes.put(columnNames[i], tokens[i]);
                }
            }
        }


        return new EQTLFeature(snp, chr, position, geneId, geneName, pValue);
    }


    /**
     * EQTL files do not have interesting headers, we are using (abusing) this method to set some initial
     * properties.  The fact we have to do this => the design needs more work.
     *
     * @param reader
     * @return
     */
    @Override
    public Object readActualHeader(LineIterator reader) {

        if (header == null) {
            header = new FeatureFileHeader();
            TrackProperties tp = new TrackProperties();
            tp.setDisplayMode(Track.DisplayMode.EXPANDED);
            header.setTrackProperties(tp);
        }
        return header;

    }

    @Override
    public boolean canDecode(String path) {

        String fn = path.toLowerCase();
        if(fn.endsWith(".gz")) fn = fn.substring(0, fn.length()-3);
        return fn.endsWith(".eqtl");
    }
}
