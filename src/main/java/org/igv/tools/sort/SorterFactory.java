package org.igv.tools.sort;

import org.igv.logging.*;
import org.igv.feature.tribble.MUTCodec;
import org.igv.gwas.GWASParser;
import org.igv.feature.gff.GFFFeatureSource;
import org.igv.util.ResourceLocator;

import java.io.File;

/**
 * Created by jrobinson on 5/10/16.
 */
public class SorterFactory {

    private static Logger log = LogManager.getLogger(SorterFactory.class);

    public static Sorter getSorter(File inputFile, File outputFile) {

        String shortFN = inputFile.getName().toLowerCase();
        if (shortFN.endsWith(".txt")) {
            shortFN = shortFN.substring(0, shortFN.length() - 4);
        }
        if (shortFN.endsWith(".cn") || shortFN.endsWith(".xcn") || shortFN.endsWith(".snp") || shortFN.endsWith(".igv")) {
            return new CNSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".sam")) {
            return new SAMSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".aligned") || shortFN.endsWith(".bed") || shortFN.endsWith(".bedgraph") || shortFN.endsWith(".bdg") ) {
            return new BedSorter(inputFile, outputFile);
        } else if (GFFFeatureSource.isGFF(shortFN)) {
            return new GFFSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".vcf")) {
            return new VCFSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".psl") || shortFN.endsWith(".pslx")) {
            return new BedSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".eqtl")) {
            return new EQTLSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".snp")) {
            return new GenericSorter(inputFile, outputFile, 1, 2);
        } else if (GWASParser.isGWASFile(shortFN)) {
            return new GWASSorter(inputFile, outputFile);
        } else if (MUTCodec.isMutationAnnotationFile(new ResourceLocator(inputFile.getAbsolutePath()))) {
            return new MUTSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".interaction")) {
            return new InteractionSorter(inputFile, outputFile);
        } else if(shortFN.endsWith(".bam")) {
            return new BAMSorter(inputFile, outputFile);
        } else if(shortFN.contains("refgene") || shortFN.contains("ncbirefseq")) {
            return new RefgeneSorter(inputFile, outputFile);
        }
        else {
            log.error("Unknown file type or sorting not supported for: " + inputFile.getName());
            return null;
        }
    }
}
