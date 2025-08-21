package org.broad.igv.htsget;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.tribble.Feature;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReader;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.VCFWrapperCodec;
import org.broad.igv.track.FeatureSource;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.*;

/**
 * A "Tribble" like feature source class for htsget variant endpoints.   Currently this only supports VCF format.  It
 * is anticipated that this class will be replaced when the htsjdk supports variant sources.
 *
 */
public class HtsgetVariantSource implements FeatureSource {

    HtsgetReader htsgetReader;
    VCFCodec codec;
    VCFHeader header;
    Map<String, String> chrAliasMap;
    int featureWindowSize;
    Genome genome;

    public HtsgetVariantSource(HtsgetUtils.Metadata metadata, Genome genome) {
        this.htsgetReader = HtsgetReader.getReader(metadata);
        this.codec = new VCFCodec();
        this.featureWindowSize = 1000;  // default
        this.chrAliasMap = new HashMap<>();
        this.genome = genome;
        init(genome);
    }

    private void init(Genome genome) {

        VCFHeader header = (VCFHeader) getHeader();
        if (genome != null) {
            List<VCFContigHeaderLine> contigsLines = header.getContigLines();
            if (contigsLines != null) {
                for (VCFContigHeaderLine l : contigsLines) {
                    String vcfChr = l.getID();
                    String genomeChr = genome.getCanonicalChrName(vcfChr);
                    chrAliasMap.put(genomeChr, vcfChr);
                }
            }
        }
    }

    @Override
    public Iterator getFeatures(String chr, int start, int end) throws IOException {

        String queryChr = chrAliasMap.containsKey(chr) ? chrAliasMap.get(chr) : chr;

        // The umccr htsget server returns bgzipped data for VCF format.  Arguably a server bug, but we
        // can handle it here.
        byte[] bytes = htsgetReader.readData(queryChr, start + 1, end);
        if (bytes != null && bytes.length >= 2 && bytes[0] == (byte) 0x1F && bytes[1] == (byte) 0x8B) {
            BlockCompressedInputStream bis = new BlockCompressedInputStream(new ByteArrayInputStream(bytes));
            bytes = bis.readAllBytes();
        }

        String data = new String(bytes);
        String[] lines = data.split("\\R");

        VCFWrapperCodec wrapperCodec = new VCFWrapperCodec(this.codec, this.genome);
        List<Feature> features = new ArrayList<>();
        for (String line : lines) {
            try {
                if (line.startsWith("#")) {
                    continue;
                } else {
                    Feature f = wrapperCodec.decode(line);
                    if (f.getEnd() < start) {
                        continue;
                    }
                    if (f.getStart() > end) {
                        break;
                    }
                    features.add(f);
                }
            } catch (Exception e) {
                // This can happen for the last feature if a partial record is returned in the query. Presumably
                // the partial record is beyond the query region.  So not actually an error.
                e.printStackTrace();
            }
        }
        return features.iterator();
    }

    @Override
    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return null;
    }

    @Override
    public int getFeatureWindowSize() {
        return featureWindowSize;
    }

    @Override
    public Object getHeader() {
        try {
            if (header == null) {

                // The UMCCR htsget server returns a bgzipped header for VCF format. Arguably a server bug,
                // but we handle it here by decompressing the header if necessary.
                byte[] bytes = htsgetReader.readHeader();
                if (bytes != null && bytes.length >= 2 && bytes[0] == (byte) 0x1F && bytes[1] == (byte) 0x8B) {
                    BlockCompressedInputStream bis = new BlockCompressedInputStream(new ByteArrayInputStream(bytes));
                    bytes = bis.readAllBytes();
                }
                String headerText = new String(bytes);

                LineIterator iter = new LineIteratorImpl(new StringLineReader(headerText));
                header = (VCFHeader) codec.readActualHeader(iter);

                // We need to parse the vcf version, its not in the header read by the codec. Perhaps an htsjdk bug?
                String formatLine = headerText.split("\\R")[0];
                VCFHeaderVersion version = VCFHeaderVersion.getHeaderVersion(formatLine);
                codec.setVCFHeader(header, version);
            }
            return header;
        } catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }
}

class StringLineReader implements LineReader {

    String[] lines;
    int lineNumber;

    public StringLineReader(String text) {
        this.lines = text.split("\\R");
        this.lineNumber = 0;
    }

    @Override
    public String readLine() throws IOException {
        if (lineNumber < lines.length) {
            return lines[lineNumber++];
        } else {
            return null;
        }

    }

    @Override
    public void close() {

    }
}
