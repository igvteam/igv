package org.igv.ucsc.bb;

import org.igv.feature.IGVFeature;
import org.igv.feature.Exon;
import org.igv.feature.PSLRecord;
import org.igv.feature.tribble.PSLCodec;
import org.igv.ucsc.bb.codecs.BBCodec;
import org.igv.ucsc.bb.codecs.BBCodecFactory;
import org.igv.ucsc.bb.codecs.BBPslCodec;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.*;

public class BBPslCodecTest {

    static final String BIGPSL_AS = "table bigPsl\n" +
            "\"bigPsl pairwise alignment\"\n" +
            "    (\n" +
            "    string chrom;       \"Reference sequence chromosome or scaffold\"\n" +
            "    uint   chromStart;  \"Start position in chromosome\"\n" +
            "    uint   chromEnd;    \"End position in chromosome\"\n" +
            "    string name;        \"Name or ID of item, ideally both human readable and unique\"\n" +
            "    uint score;         \"Score (0-1000)\"\n" +
            "    char[1] strand;     \"+ or - indicates whether the query aligns to the + or - strand on the reference\"\n" +
            "    uint thickStart;    \"Start of where display should be thick (start codon)\"\n" +
            "    uint thickEnd;      \"End of where display should be thick (stop codon)\"\n" +
            "    uint reserved;      \"RGB value (use R,G,B string in input file)\"\n" +
            "    int blockCount;     \"Number of blocks\"\n" +
            "    int[blockCount] blockSizes; \"Comma separated list of block sizes\"\n" +
            "    int[blockCount] chromStarts; \"Start positions relative to chromStart\"\n" +
            "    uint    oChromStart;\"Start position in other chromosome\"\n" +
            "    uint    oChromEnd;  \"End position in other chromosome\"\n" +
            "    char[1] oStrand;    \"+ or -, - means that psl was reversed into BED-compatible coordinates\"\n" +
            "    uint    oChromSize; \"Size of other chromosome.\"\n" +
            "    int[blockCount] oChromStarts; \"Start positions relative to oChromStart or from oChromStart+oChromSize depending on strand\"\n" +
            "    lstring  oSequence;  \"Sequence on other chrom (or edit list, or empty)\"\n" +
            "    string   oCDS;       \"CDS in NCBI format\"\n" +
            "    uint    chromSize;  \"Size of target chromosome\"\n" +
            "    uint match;        \"Number of bases matched.\"\n" +
            "    uint misMatch; \" Number of bases that don't match \"\n" +
            "    uint repMatch; \" Number of bases that match but are part of repeats \"\n" +
            "    uint nCount;   \" Number of 'N' bases \"\n" +
            "    uint seqType;    \"0=empty, 1=nucleotide, 2=amino_acid\"\n" +
            "    )\n";

    // Second record of test/data/psl/humanBlat.psl
    static final String PSL_LINE = "27\t0\t0\t0\t0\t0\t1\t145\t+\tYourSeq\t119\t92\t119\tchr6\t170805979\t84433944\t84434116\t2\t18,9,\t92,110,\t84433944,84434107,";

    @Test
    public void testDecode() {

        BBCodec codec = BBCodecFactory.getCodec(BIGPSL_AS, 12);
        assertTrue(codec instanceof BBPslCodec);

        // Same alignment as PSL_LINE in bigPsl form (fields 4-25)
        String rest = "YourSeq\t1000\t+\t84433944\t84434116\t0\t2\t18,9,\t0,163,\t92\t119\t+\t119\t92,110,\t\t\t170805979\t27\t0\t0\t0\t0";
        IGVFeature f = codec.decode(new BedData("chr6", 84433944, 84434116, rest));

        assertTrue(f instanceof PSLRecord);
        PSLRecord record = (PSLRecord) f;
        assertEquals("chr6", record.getChr());
        assertEquals(84433944, record.getStart());
        assertEquals(84434116, record.getEnd());
        assertEquals("YourSeq", record.getName());
        assertEquals(27, record.getMatch());
        assertEquals(119, record.getqSize());
        assertEquals(170805979, record.getTSize());
        assertEquals(0, record.getQGapCount());
        assertEquals(0, record.getQGapBases());
        assertEquals(1, record.getTGapCount());
        assertEquals(145, record.getTGapBases());

        List<Exon> exons = record.getExons();
        assertEquals(2, exons.size());
        assertEquals(84433944, exons.get(0).getStart());
        assertEquals(18, exons.get(0).getLength());
        assertEquals(84434107, exons.get(1).getStart());
        assertEquals(9, exons.get(1).getLength());

        // Identical to the ascii record
        PSLRecord ascii = new PSLCodec().decode(PSL_LINE);
        assertEquals(ascii.getScore(), record.getScore(), 0.01);
        assertEquals(ascii.getDescription(), record.getDescription());
        assertEquals(ascii.getText(), record.getText());
    }
}
