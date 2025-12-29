package org.igv.bedpe;

import org.igv.Globals;
import org.igv.feature.genome.Genome;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.util.ParsingUtils;
import org.igv.util.ResourceLocator;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.List;
import java.util.*;

/**
 * Created by jrobinso on 6/29/18.
 */
public class InteractParser {

    private static Logger log = LogManager.getLogger(InteractParser.class);

    public static List<BedPE> parse(ResourceLocator locator, Genome genome) throws IOException {

        List<BedPE> features = new ArrayList<>();
        BufferedReader br = null;

        try {
            br = ParsingUtils.openBufferedReader(locator.getPath());
            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                if (nextLine.startsWith("#") || nextLine.startsWith("track") || nextLine.startsWith("browser")) continue;
                String [] tokens = Globals.whitespacePattern.split(nextLine);
                InteractFeature f = InteractFeature.fromTokens(tokens, genome);
                if(f != null) {
                    features.add(f);
                }
            }
            return  features;
        } finally {
            br.close();
        }
    }

}
