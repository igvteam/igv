package org.broad.igv.feature.genome.load;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.Sequence;
import org.broad.igv.feature.genome.SequenceWrapper;
import org.broad.igv.feature.genome.fasta.FastaBlockCompressedSequence;
import org.broad.igv.feature.genome.fasta.FastaIndexedSequence;
import org.broad.igv.feature.genome.fasta.FastaUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.Utilities;

import java.io.File;
import java.io.IOException;

public class FastaGenomeLoader extends GenomeLoader {

    private String genomePath;


    public FastaGenomeLoader(String genomePath) {
        this.genomePath = genomePath;
    }

    @Override
    public Genome loadGenome() throws IOException {

        String fastaPath = null;
        String fastaIndexPath = null;
        if (genomePath.endsWith(".fai")) {
            fastaPath = genomePath.substring(0, genomePath.length() - 4);
            fastaIndexPath = genomePath;
        } else {
            fastaPath = genomePath;
            fastaIndexPath = genomePath + ".fai";
        }

        if (!FileUtils.resourceExists(fastaIndexPath)) {
            //Have to make sure we have a local copy of the fasta file
            //to index it
            if (!FileUtils.isRemote(fastaPath)) {
                fastaIndexPath = fastaPath + ".fai";
                FastaUtils.createIndexFile(fastaPath, fastaIndexPath);
            }
        }

        String id = fastaPath;
        String name;
        if (HttpUtils.isRemoteURL(fastaPath)) {
            name = Utilities.getFileNameFromURL(fastaPath);
        } else {
            File file = new File(fastaPath);
            if (!file.exists()) {
                throw new IOException(fastaPath + " does not exist, could not load genome");
            }
            name = file.getName();
        }

        FastaIndexedSequence fastaSequence = fastaPath.endsWith(".gz") ?
                new FastaBlockCompressedSequence(fastaPath) :
                new FastaIndexedSequence(fastaPath);
        Sequence sequence = new SequenceWrapper(fastaSequence);
        return new Genome(id, name, sequence, true);
    }


}
