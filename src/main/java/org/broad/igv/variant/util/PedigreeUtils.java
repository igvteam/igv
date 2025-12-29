package org.broad.igv.variant.util;

import org.broad.igv.track.AttributeManager;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.TrackGroupEvent;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;

/**
 * @author jrobinso
 * @date Apr 25, 2011
 */
public class PedigreeUtils {


    /**
     * Parses a pedigree file containing trios of child, father, mother.
     *
     * @param path
     * @throws IOException
     */

    public static void parseTrioFile(String path) throws IOException {

        BufferedReader reader = null;


        try {
            reader = ParsingUtils.openBufferedReader(path);
            String nextLine;
            int familyNumber = 1;
            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = nextLine.split("\\s+");
                if (tokens.length == 3) {
                    String family = String.valueOf(familyNumber);
                    for (String member : tokens) {
                        AttributeManager.getInstance().addAttribute(member, "Family", family);
                    }
                    String child = tokens[0];
                    String father = tokens[1];
                    String mother = tokens[2];
                    AttributeManager.getInstance().addAttribute(child, "Member", "child");
                    AttributeManager.getInstance().addAttribute(child, "Relation", "child");
                    AttributeManager.getInstance().addAttribute(father, "Member", "father");
                    AttributeManager.getInstance().addAttribute(father, "Relation", "parent");
                    AttributeManager.getInstance().addAttribute(mother, "Member", "mother");
                    AttributeManager.getInstance().addAttribute(mother, "Relation", "parent");
                    familyNumber++;
                }
            }

        } finally {
            if (reader != null) reader.close();
        }

    }
}
