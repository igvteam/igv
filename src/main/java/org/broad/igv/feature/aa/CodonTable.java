package org.broad.igv.feature.aa;

import com.google.common.base.Objects;
import com.google.common.collect.Table;
import com.google.common.collect.TreeBasedTable;
import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonParser;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;

/**
 * Store information about current codon translation table.
 * Intended to be loaded from external resource, and then never modified.
 * To that end, collections contained here are set to be unmodifiable
 */
public class CodonTable {


    private static final String[] BASE_SEQUENCES = {"TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
            "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
            "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"};


    private final Integer id;
    private final List<String> names;
    private final Set<AminoAcid> starts;
    private final Map<String, AminoAcid> codonMap;
    private Set<String> altStartCodons;



    /**
     * Get the amino acid represented by this codon
     *
     * @param codon
     * @return
     */
    public AminoAcid getAminoAcid(String codon) {
        if (codon.length() != 3) {
            throw new IllegalArgumentException("Codon must be length 3: " + codon);
        }

        AminoAcid aa = codonMap.get(codon);
        if (aa == null) {
            return AminoAcid.NULL_AMINO_ACID;
        }
        return aa;
    }

     CodonTable(String path, int id, List<String> names, Set<AminoAcid> starts, Map<String, AminoAcid> codonMap) {
        this.id = id;
        this.names = Collections.unmodifiableList(names);
        this.starts = Collections.unmodifiableSet(starts);
        this.codonMap = Collections.unmodifiableMap(codonMap);
    }

     static CodonTable createFromJSON(String sourcePath, JsonObject jsonObject) throws JsonParseException {
        int id = jsonObject.get("id").getAsInt();

        JsonArray jsonnames = jsonObject.get("name").getAsJsonArray();
        List<String> names = new ArrayList<String>(jsonnames.size());
        for (int nn = 0; nn < jsonnames.size(); nn++) {
            names.add(jsonnames.get(nn).getAsString());
        }

        //Data is written as several long strings which line up
        String aas = jsonObject.get("ncbieaa").getAsString();
        String startString = jsonObject.get("sncbieaa").getAsString();

        CodonTable codonTable = build(sourcePath, id, names, aas, startString);
        if (jsonObject.has("altStartCodons")) {
            JsonArray a = jsonObject.get("altStartCodons").getAsJsonArray();
            Set<String> altStartCodons = new HashSet<>();
            for (int i = 0; i < a.size(); i++) {
                altStartCodons.add(a.get(i).getAsString());
            }
            codonTable.altStartCodons = altStartCodons;
        }

        return codonTable;
    }

    private static CodonTable build(String sourcePath, int id, List<String> names, String aas, String startString) {

        String base1 = BASE_SEQUENCES[0];
        String base2 = BASE_SEQUENCES[1];
        String base3 = BASE_SEQUENCES[2];

        checkLengths(base1, base2, base3, aas, startString);


        Map<String, AminoAcid> codonMap = new HashMap<String, AminoAcid>(aas.length());
        Set<AminoAcid> starts = new HashSet<AminoAcid>(aas.length());

        for (int cc = 0; cc < aas.length(); cc++) {
            String codon = base1.substring(cc, cc + 1) + base2.substring(cc, cc + 1) + base3.substring(cc, cc + 1);
            AminoAcid aa = AminoAcidManager.AANameMap.get(aas.substring(cc, cc + 1));

            codonMap.put(codon, aa);

            if (startString.charAt(cc) == 'M') {
                starts.add(aa);
            }
        }

        return new CodonTable(sourcePath, id, names, starts, codonMap);
    }

    private static void checkLengths(String... values) {
        int length = values[0].length();
        assert length == 64;
        for (int v = 1; v < values.length; v++) {
            if (values[v].length() != length) {
                String msg = "Amino acid and codon strings must all be the same length.";
                msg += "Expected length " + length + ", found length " + values[v].length();
                throw new InputMismatchException(msg);
            }
        }
    }

    public int getId() {
        return id;
    }

    public String getDisplayName() {
        return names.get(0);
    }

    public Set<AminoAcid> getStarts() {
        return starts;
    }

    Map<String, AminoAcid> getCodonMap() {
        return codonMap;
    }

    public Set<String> getAltStartCodons() {
        return altStartCodons;
    }

    @Override
    public boolean equals(Object object) {
        if (object instanceof CodonTable) {
            CodonTable other = (CodonTable) object;
            return com.google.common.base.Objects.equal(this.id, other.id) &&
                    com.google.common.base.Objects.equal(this.names, other.names) &&
                    com.google.common.base.Objects.equal(this.starts, other.starts) &&
                    com.google.common.base.Objects.equal(this.codonMap, other.codonMap);
        }
        return false;
    }

    @Override
    public int hashCode() {
        return Objects.hashCode(this.id, this.names, this.starts, this.codonMap);
    }

}
