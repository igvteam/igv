package org.igv.data;

import org.igv.track.TrackType;
import org.igv.util.ParsingUtils;

/**
 * @author jrobinso
 * @date Apr 23, 2011
 */
public class AbstractDataParser {

      protected void parseDirective(String comment, IGVDataset dataset) {

        String tmp = comment.substring(1, comment.length());
        if (tmp.startsWith("track")) {
            ParsingUtils.parseTrackLine(tmp, dataset.getTrackProperties());

        } else {
            String[] tokens = tmp.split("=");
            if (tokens.length != 2) {
                return;
            }

            String key = tokens[0].trim().toLowerCase();
            if (key.equals("name")) {
                dataset.setName(tokens[1].trim());
            } else if (key.equals("type")) {

                try {
                    dataset.setTrackType(TrackType.valueOf(tokens[1].trim().toUpperCase()));
                } catch (Exception exception) {

                    // Ignore
                }
            }
        }
    }
}
