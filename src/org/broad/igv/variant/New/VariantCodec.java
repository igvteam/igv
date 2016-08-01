package org.broad.igv.variant.New;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.StringUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Created by jrobinso on 7/29/16.
 */
public class VariantCodec extends AsciiFeatureCodec<Variant> {

    private static Logger log = Logger.getLogger(VariantCodec.class);

    Genome genome;
    VCFHeader header;

    protected VariantCodec(Genome genome, ResourceLocator locator) throws IOException {
        super(Variant.class);
        this.genome = genome;

        BufferedReader reader = ParsingUtils.openBufferedReader(locator);
        header = parseHeader(reader);
    }

    @Override
    public Variant decode(String line) {


        String[] tokens = Globals.tabPattern.split(line);

        if (tokens.length > 7) {

            Variant variant = new Variant();
            variant.chr = genome == null ? tokens[0] : genome.getCanonicalChrName(tokens[0]);
            variant.pos = Integer.parseInt(tokens[1]) - 1;
            variant.names = tokens[2];    // id in VCF
            variant.referenceBases = tokens[3];
            variant.alternateBases = tokens[4];
            variant.quality = Integer.parseInt(tokens[5]);
            variant.filter = tokens[6];
            variant.info = tokens[7];      // TODO -- parse special reserved values, e.g. AF

            String[] altTokens = variant.alternateBases.split(",");
            variant.alleles = new String[altTokens.length + 1];
            variant.alleles[0] = variant.referenceBases;
            for (int i = 0; i < altTokens.length; i++) {
                variant.alleles[i + 1] = altTokens[i];
            }

            computeStart(variant);

            // TODO -- Calls

            return variant;

        } else {
            // invalid line
            return null;
        }
    }

    private static void computeStart(Variant variant) {

        if (variant.alleles.length > 1) {

            variant.start = variant.pos;
            variant.end = 0;

            int refBasesLength = variant.alleles[0].length();

            for (int i = 1; i < variant.alleles.length; i++) {

                int altLength = variant.alleles[i].length();

                if (altLength > 0) {

                    int s, e;
                    int diff = refBasesLength - altLength;

                    if (diff > 0) {
                        // deletion, assume left padded
                        s = variant.pos + altLength;
                        e = s + diff;
                    } else if (diff < 0) {
                        // Insertion, assume left padded, insertion begins to "right" of last ref base
                        s = variant.pos + refBasesLength;
                        e = s + 1;     // Insertion between s & 3
                    } else {
                        // Substitution, SNP if seq.length == 1
                        s = variant.pos;
                        e = s + altLength;
                    }
                    // variant.alleles.push({allele: alt, start: s, end: e});
                    variant.start = Math.min(variant.start, s);
                    variant.end = Math.max(variant.end, e);
                }

            }
        } else {
            // Is this even legal VCF?  (NO alt alleles)
            variant.start = variant.pos - 1;
            variant.end = variant.pos;
        }
    }


    @Override
    public Object readActualHeader(LineIterator lineIterator) {
        return null;   // Who actuall calls this?
    }

    @Override
    public boolean canDecode(String s) {
        return false;
    }


    static class VCFHeader {

        String version;
        Map<String, Map<String, String>> infoFields;
        Map<String, Map<String, String>> formatFields;
        Map<String, Map<String, String>> filterFields;
        Map<Integer, String> callSetNames;

        public VCFHeader() {
            infoFields = new HashMap<>();
            formatFields = new HashMap<>();
            filterFields = new HashMap<>();
        }

        public void addMetaInfo(String type, String id, Map<String, String> info) {

            Map<String, Map<String, String>> dict;
            if ("##INFO".equals(type)) {
                dict = infoFields;
            } else if ("##FORMAT".equals(type)) {
                dict = formatFields;
            } else if ("##FILTER".equals(type)) {
                dict = filterFields;
            } else {
                // ignored
                return;
            }

            dict.put(type, info);

        }


    }

    private VCFHeader parseHeader(BufferedReader reader) throws IOException {

        VCFHeader header = new VCFHeader();

        String line = reader.readLine();
        // First line must be file format
        if (line != null && line.startsWith("##fileformat")) {
            header.version = line.substring(13);
        } else {
            throw new Error("Invalid VCF file: missing fileformat line");
        }

        while ((line = reader.readLine()) != null) {

            if (line.startsWith("#")) {

                if (line.startsWith("##")) {

                    if (line.startsWith("##INFO") || line.startsWith("##FILTER") || line.startsWith("##FORMAT")) {

                        int ltIdx = line.indexOf("<");
                        int gtIdx = line.lastIndexOf(">");

                        if (!(ltIdx > 2 && gtIdx > 0)) {
                            log.error("Malformed VCF header line: " + line);
                            continue;
                        }

                        //##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency based on Flow Evaluator observation counts">
                        // ##FILTER=<ID=NOCALL,Description="Generic filter. Filtering details stored in FR info tag.">
                        // ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele frequency based on Flow Evaluator observation counts">

                        String type = line.substring(2, ltIdx - 1);

                        List<String> tokens = StringUtils.breakQuotedString(line.substring(ltIdx + 1, gtIdx - 1), ',');

                        String id = null;
                        Map<String, String> values = new HashMap<>();
                        for (String token : tokens) {

                            String[] kv = token.split("=");
                            if (kv.length > 1) {
                                if ("ID".equals(kv[0])) {
                                    id = kv[1];
                                } else {
                                    values.put(kv[0], kv[1]);
                                }
                            }
                        }
                        ;

                        if (id != null) {
                            header.addMetaInfo(type, id, values);
                        }
                    } else {
                        // Ignoring other ## header lines
                    }
                } else if (line.startsWith("#CHROM")) {

                    String[] tokens = line.split("\t");

                    if (tokens.length > 8) {

                        // call set names -- create column # -> cs name table
                        header.callSetNames = new HashMap<>();
                        for (int j = 9; j < tokens.length; j++) {
                            header.callSetNames.put(j, tokens[j]);
                        }
                    }
                    break;
                }

            } else {
                break;
            }

        }

        return header;

    }
//
//    function extractCallFields(tokens) {
//
//        var callFields = {
//                genotypeIndex:-1,
//                genotypeLikelihoodIndex:-1,
//                phasesetIndex:-1,
//                fields:tokens
//        },
//        i;
//
//        for (i = 0; i < tokens.length; i++) {
//
//            if ("GT" == = tokens[i]) {
//                callFields.genotypeIndex = i;
//            } else if ("GL" == = tokens[i]) {
//                callFields.genotypeLikelihoodIndex = i;
//            } else if ("PS" == = tokens[i]) {
//                callFields.phasesetIndex = i;
//            }
//        }
//        return callFields;
//
//    }
}
