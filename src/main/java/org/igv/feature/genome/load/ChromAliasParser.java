package org.igv.feature.genome.load;

import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.ui.util.MessageUtils;
import org.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class ChromAliasParser {

    private static Logger log = LogManager.getLogger(ChromAliasParser.class);

    /**
     * Load user-defined chromosome aliases.
     *
     * @param path
     * @return
     */
    public static List<List<String>> loadChrAliases(String path) {
        try (BufferedReader br = ParsingUtils.openBufferedReader(path)) {
            return loadChrAliases(br);
        } catch (IOException e) {
            log.error("Error loading chr alias table", e);
            MessageUtils.showMessage("<html>Error loading chromosome alias table.  Aliases will not be available<br>" +
                    e.toString());
            return null;
        }
    }

    public static List<List<String>> loadChrAliases(BufferedReader br) throws IOException {
        String nextLine = "";
        List<List<String>> synonymList = new ArrayList<List<String>>();
        while ((nextLine = br.readLine()) != null) {
            String[] tokens = nextLine.split("\t");
            if (tokens.length > 1) {
                List<String> synonyms = new ArrayList<String>();
                for (String t : tokens) {
                    String syn = t.trim();
                    if (t.length() > 0) synonyms.add(syn.trim());
                }
                synonymList.add(synonyms);
            }
        }
        return synonymList;
    }
}
