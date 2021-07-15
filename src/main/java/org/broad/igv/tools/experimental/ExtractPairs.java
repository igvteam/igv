package org.broad.igv.tools.experimental;

import htsjdk.samtools.*;
import org.broad.igv.sam.cram.IGVReferenceSource;
import org.broad.igv.util.HttpUtils;

import java.io.*;
import java.net.URL;

public class ExtractPairs {

    private static final int READ_PAIRED_FLAG = 0x1;
    private static final int PROPER_PAIR_FLAG = 0x2;
    private static final int READ_UNMAPPED_FLAG = 0x4;
    private static final int MATE_UNMAPPED_FLAG = 0x8;
    private static final int READ_STRAND_FLAG = 0x10;
    protected static final int MATE_STRAND_FLAG = 0x20;
    private static final int FIRST_OF_PAIR_FLAG = 0x40;
    private static final int SECOND_OF_PAIR_FLAG = 0x80;
    private static final int NOT_PRIMARY_ALIGNMENT_FLAG = 0x100;
    private static final int READ_FAILS_VENDOR_QUALITY_CHECK_FLAG = 0x200;
    private static final int DUPLICATE_READ_FLAG = 0x400;
    private static final int SUPPLEMENTARY_ALIGNMENT_FLAG = 0x800;

    public static void main(String[] args) throws IOException {

        String inputFile = "/Volumes/GoogleDrive/Shared drives/IGV/MB275/MB275.discordant.namesorted.bam";
        String outputFile = "/Volumes/GoogleDrive/Shared drives/IGV/MB275/pairs2.txt";
        extract(inputFile, outputFile);
    }

    // * str1 chr1 pos1 frag1 str2 chr2 pos2 frag2

    static void extract(String alignmentFile, String outputFile) throws IOException {

        SamReader reader = null;
        SAMRecordIterator iter = null;

        PrintWriter pw = null;

        try {
            pw = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));  //new PrintWriter(new OutputStreamWriter(System.out)); //
            reader = getSamReader(alignmentFile);
            iter = reader.iterator();

            SAMRecord lastAlignment = null;

            while (iter != null && iter.hasNext()) {
                SAMRecord alignment = iter.next();
                if (filter(alignment)) {
                    if (lastAlignment == null || !(alignment.getReadName().equals(lastAlignment.getReadName()))) {
                        lastAlignment = alignment;
                    } else {
                        // pair
                        // Not looking for orientation discordancy, filter
                        if (!lastAlignment.getContig().equals(alignment.getContig()) ||
                                Math.abs(lastAlignment.getStart() - alignment.getStart()) > 1000) {

                            SAMRecord first = lastAlignment.getContig().compareTo(alignment.getContig()) < 0 ? lastAlignment : alignment;
                            SAMRecord second = lastAlignment.getContig().compareTo(alignment.getContig()) < 0 ? alignment : lastAlignment;
                            pw.println(posString(first) + "\t0\t" + posString(second) + "\t0");
                        }
                        lastAlignment = null;
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            pw.flush();
            pw.close();
            iter.close();
            reader.close();
        }
    }

    static boolean filter(SAMRecord alignment) {
        return alignment.getMappingQuality() > 0;
    }

    static String posString(SAMRecord a) {
        return ((a.getFlags() & READ_STRAND_FLAG) == 0 ? "0" : "1") + "\t" +
                a.getContig() + "\t" + ((a.getStart() + a.getEnd()) / 2);
    }

    private static SamReader getSamReader(String path) throws IOException {

        boolean isLocal = !(path.startsWith("http"));
        final SamReaderFactory factory = SamReaderFactory.makeDefault().
                referenceSource(new IGVReferenceSource()).
                validationStringency(ValidationStringency.SILENT);
        SamInputResource resource;

        if (isLocal) {
            resource = SamInputResource.of(new File(path));
        } else {
            URL url = new URL(path);
            resource = SamInputResource.of(new BufferedInputStream(HttpUtils.getInstance().openConnectionStream(url)));
        }

        return factory.open(resource);

    }

}
