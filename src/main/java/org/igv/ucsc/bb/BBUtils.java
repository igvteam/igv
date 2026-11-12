package org.igv.ucsc.bb;

import org.igv.Globals;

import java.util.ArrayList;
import java.util.List;

public class BBUtils {

    public static ASTable parseAutosql(String str) {

        List<ASField> fields = new ArrayList<>();
        boolean startDecoding = false;
        String table = "";
        String[] lines = str.trim().split("\\R");
        for (String line : lines) {
            if (line.startsWith("#")) {
                continue;
            }
            line = line.trim();
            if (line.length() > 0) {
                if (line.startsWith("table")) {
                    table = line.split("\\s+")[1].trim();
                } else if (line.startsWith("(")) {
                    startDecoding = true;
                } else if (line.startsWith(")")) {
                    break;
                } else if (startDecoding) {
                    int idx = line.indexOf(";");
                    if (idx > 0) {
                        String[] tokens = Globals.whitespacePattern.split(line.substring(0, idx));
                        String description = line.substring(idx + 1).replace("\"", "").trim();
                        fields.add(new ASField(tokens[0], tokens[1], description));
                    }
                }
            }
        }

        return new ASTable(table, fields);
    }

    private static final String[][] BED_FIELD_NAMES = {{"chrom"}, {"chromStart"}, {"chromEnd"}, {"name"}, {"score"},
            {"strand"}, {"thickStart"}, {"thickEnd"}, {"reserved", "itemRgb", "color"}, {"blockCount"},
            {"blockSizes"}, {"chromStarts", "blockStarts"}};

    /**
     * Infer the number of standard BED fields for files that do not record it (definedFieldCount == 0 in the header,
     * as written by early versions of bedToBigBed).  With autosql, count the leading fields with standard BED names.
     * Without autosql all fields are assumed to be standard.
     */
    public static int inferDefinedFieldCount(String autosql, int fieldCount) {
        if (autosql == null || autosql.length() == 0) {
            return Math.max(3, Math.min(fieldCount, BED_FIELD_NAMES.length));
        }
        List<ASField> fields = parseAutosql(autosql).fields;
        int n = 3;
        while (n < fields.size() && n < BED_FIELD_NAMES.length && n < fieldCount) {
            String name = fields.get(n).name;
            boolean match = false;
            for (String bedName : BED_FIELD_NAMES[n]) {
                if (bedName.equalsIgnoreCase(name)) {
                    match = true;
                    break;
                }
            }
            if (!match) break;
            n++;
        }
        return n;
    }

    public static class ASTable {
        public String name;
        public List<ASField> fields;

        public ASTable(String name, List<ASField> fields) {
            this.name = name;
            this.fields = fields;
        }
    }

    public static class ASField {
        public String type;
        public String name;
        public String description;

        public ASField(String type, String name, String description) {
            this.type = type;
            this.name = name;
            this.description = description;
        }
    }

}


