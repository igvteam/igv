package org.igv.ui.genome;

import org.igv.feature.genome.GenomeDownloadUtils;
import org.igv.feature.genome.load.GenomeConfig;
import org.igv.feature.genome.load.JsonGenomeLoader;

import java.io.File;
import java.io.IOException;

public class GenomeUtils {

    public static GenomeListItem updateGenome(GenomeListItem hostedItem) throws IOException {

        String id = hostedItem.getId();
        String path = hostedItem.getPath();

        if (path.endsWith(".json")) {
            JsonGenomeLoader loader = new JsonGenomeLoader(hostedItem.getPath());
            GenomeConfig config = loader.loadGenomeConfig();
            File localFile = GenomeDownloadUtils.downloadGenome(config, false, false);

            if (localFile != null) {
                return new GenomeListItem(hostedItem.getDisplayableName(), localFile.getAbsolutePath(), id);
            } else {
                return new GenomeListItem(hostedItem.getDisplayableName(), hostedItem.getPath(), id);
            }
        }

        return hostedItem;
    }

    public static boolean isDeprecated(GenomeConfig config) {
        // TwoBit URLs are not deprecated and take precedence over fastaURL
        if (config.twoBitURL != null && !config.twoBitURL.isEmpty()) {
            return false;
        }

        // Check if fastaURL is deprecated
        String fastaURL = config.fastaURL;
        return isDeprecated(fastaURL);
    }

    public static boolean isDeprecated(String fastaURL) {
        return fastaURL.contains("amazonaws") &&
                (fastaURL.contains("igv.org.genomes") ||
                        fastaURL.contains("igv.broadinstitute.org") ||
                        fastaURL.contains("igv-genepattern-org"));
    }
}
