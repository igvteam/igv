package org.igv.feature.genome.load;

import org.igv.Globals;
import org.igv.feature.Chromosome;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 4/16/13
 *         Time: 2:25 PM
 */
public class ChromSizesParser {

    private static Logger log = LogManager.getLogger(ChromSizesParser.class);

    public static List<Chromosome> parse(String path) throws IOException {

        BufferedReader br = null;

        try {
            br = ParsingUtils.openBufferedReader(path);
            List<Chromosome> chromosomes = new ArrayList<Chromosome>();
            String nextLine;
            int idx = 0;
            while ((nextLine = br.readLine()) != null) {

                String[] tokens = Globals.whitespacePattern.split(nextLine);
                if (tokens.length >= 2) {
                    String chr = tokens[0];
                    int length = Integer.parseInt(tokens[1]);
                    chromosomes.add(new Chromosome(idx, chr, length));
                    idx++;
                } else {
                    log.warn("Unexpected # of tokens at line: " + nextLine);
                }

            }
            return chromosomes;
        } finally {
            if (br != null) br.close();
        }

    }


}
