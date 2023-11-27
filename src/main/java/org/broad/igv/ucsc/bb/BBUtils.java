package org.broad.igv.ucsc.bb;

import org.broad.igv.Globals;

import java.util.ArrayList;
import java.util.List;

public class BBUtils {

    public static ASTable parseAutosql(String str) {

        List<ASField> fields = new ArrayList<>();
        boolean startDecoding = false;
        String table = "";
        String[] lines = str.trim().split("\\R");
        for (String line : lines) {
            line = line.trim();
            if (line.startsWith("table")) {
                table = line.split("\\s+")[1].trim();
            } else if (line.startsWith("(")) {
                startDecoding = true;
            } else if (line.startsWith(")")) {
                break;
            } else if (startDecoding) {
                if (line.length() > 0) {
                    //                "    string chrom;       \"Reference sequence chromosome or scaffold\"\n" +
                    int idx = line.indexOf(";");
                    String[] tokens = Globals.whitespacePattern.split(line.substring(0, idx));
                    String description = line.substring(idx + 1).replace("\"", "").trim();
                    fields.add(new ASField(tokens[0], tokens[1], description));
                }
            }
        }

        return new ASTable(table, fields);
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


