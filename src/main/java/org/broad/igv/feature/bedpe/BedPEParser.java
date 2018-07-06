package org.broad.igv.feature.bedpe;

import org.apache.log4j.Logger;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by jrobinso on 6/29/18.
 */
public class BedPEParser {

    private static Logger log = Logger.getLogger(BedPEParser.class);

    public static List<BedPEFeature> parse(String file) throws IOException {

        List<BedPEFeature> features = new ArrayList<>();

        BufferedReader br = null;

        br = ParsingUtils.openBufferedReader(file);

        String nextLine;
        while ((nextLine = br.readLine()) != null) {

            String[] tokens = ParsingUtils.WHITESPACE_PATTERN.split(nextLine);

            if (tokens.length < 6) {
                log.info("Skipping line: " + nextLine);
                continue;
            }

            BedPEFeature feature = new BedPEFeature();
            feature.chr1 = tokens[0];
            feature.start1 = Integer.parseInt(tokens[1]);
            feature.end1 = Integer.parseInt(tokens[2]);
            feature.chr2 = tokens[3];
            feature.start2 = Integer.parseInt(tokens[4]);
            feature.end2 = Integer.parseInt(tokens[5]);

            if (tokens.length > 6) {
                feature.name = tokens[6];
            }

            if (tokens.length > 7) {
                feature.score = tokens[7];
            }

            // Skipping remaining fields for now

            features.add(feature);

        }

        return features;


    }


}
