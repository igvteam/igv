package org.broad.igv.htsget;

import org.apache.log4j.Logger;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.htsget.HtsgetClass;
import htsjdk.samtools.util.htsget.HtsgetFormat;
import htsjdk.samtools.util.htsget.HtsgetRequest;

import java.net.URI;
import java.net.URISyntaxException;


public class Htsget {
    private static Logger log = Logger.getLogger(Htsget.class);

    private static final String endpoint = "https://htsget.ga4gh.org/reads/10X_P4_0_possorted_genome.bam";

    public void htsgetRequest() {
        try {
            final HtsgetRequest req = new HtsgetRequest(new URI(endpoint + 1))
                    .withFormat(HtsgetFormat.BAM)
                    .withDataClass(HtsgetClass.body)
                    .withInterval(new Interval("chr1", 1, 16));
            final String query = req.toURI().toString();
        } catch (URISyntaxException e) {
            log.error(e.getReason());
        }
    }
}