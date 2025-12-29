package org.igv.feature.genome.load;

import org.igv.feature.genome.Genome;
import org.igv.feature.genome.fasta.FastaUtils;
import org.igv.util.FileUtils;
import org.igv.util.HttpUtils;
import org.igv.util.Utilities;

import java.io.File;
import java.io.IOException;

/**
 * Loader for a minimal genome of only reference sequence, which can be defined by a fasta or 2bit file.
 */
public class FastaGenomeLoader extends GenomeLoader {

    private String genomePath;


    public FastaGenomeLoader(String genomePath) {
        this.genomePath = genomePath;
    }

    @Override
    public Genome loadGenome() throws IOException {

        GenomeConfig config = new GenomeConfig();
        String id = genomePath;
        String name;
        if (HttpUtils.isRemoteURL(genomePath)) {
            name = Utilities.getFileNameFromURL(genomePath);
        } else {
            File file = new File(genomePath);
            if (!file.exists()) {
                throw new IOException(genomePath + " does not exist, could not load genome");
            }
            name = file.getName();
        }


        config.id = (id);
        config.setName(name);

        if (genomePath.endsWith(".2bit")) {
            config.twoBitURL =(genomePath);

        } else {
            String fastaPath;
            String fastaIndexPath;
            if (genomePath.endsWith(".fai")) {
                fastaPath = genomePath.substring(0, genomePath.length() - 4);
                fastaIndexPath = genomePath;
            } else {
                fastaPath = genomePath;
                fastaIndexPath = genomePath + ".fai";
            }

            // If fasta is a local file create index, if needed
            if (!FileUtils.isRemote(fastaPath) && !FileUtils.resourceExists(fastaIndexPath)) {
                FastaUtils.createIndexFile(fastaPath, fastaIndexPath);
            }

            config.fastaURL = (fastaPath);
            config.indexURL = (fastaIndexPath);
        }

        return new Genome(config);
    }


}
