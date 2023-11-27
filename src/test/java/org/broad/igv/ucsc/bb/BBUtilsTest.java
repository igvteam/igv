/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.ucsc.bb;

import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;

/**
 * @author jrobinso
 * @date Jun 19, 2011
 */
public class BBUtilsTest {


    @Test
    public void testAutosql()  {

        String autosql = "table bigCat\n" +
                "\"bigCat gene models\"\n" +
                "    (\n" +
                "    string chrom;       \"Reference sequence chromosome or scaffold\"\n" +
                "    uint   chromStart;  \"Start position in chromosome\"\n" +
                "    uint   chromEnd;    \"End position in chromosome\"\n" +
                "    string name;        \"Name\"\n" +
                "    uint score;         \"Score (0-1000)\"\n" +
                "    char[1] strand;     \"+ or - for strand\"\n" +
                "    uint thickStart;    \"Start of where display should be thick (start codon)\"\n" +
                "    uint thickEnd;      \"End of where display should be thick (stop codon)\"\n" +
                "    uint reserved;       \"RGB value (use R,G,B string in input file)\"\n" +
                "    int blockCount;     \"Number of blocks\"\n" +
                "    int[blockCount] blockSizes; \"Comma separated list of block sizes\"\n" +
                "    int[blockCount] chromStarts; \"Start positions relative to chromStart\"\n" +
                "    string name2;       \"Gene name\"\n" +
                "    string cdsStartStat; \"Status of CDS start annotation\"\n" +
                "    string cdsEndStat;   \"Status of CDS end annotation\"\n" +
                "    int[blockCount] exonFrames; \"Exon frame {0,1,2}, or -1 if no frame for exon\"\n" +
                "    string txId; \"Transcript ID\"\n" +
                "    string type;        \"Transcript type\"\n" +
                "    string geneName;    \"Gene ID\"\n" +
                "    string geneType;    \"Gene type\"\n" +
                "    string sourceGene;    \"Source gene ID\"\n" +
                "    string sourceTranscript;    \"Source transcript ID\"\n" +
                "    string alignmentId;  \"Alignment ID\"\n" +
                "    lstring alternativeSourceTranscripts;    \"Alternative source transcripts\"\n" +
                "    lstring Paralogy;    \"Paralogous alignment IDs\"\n" +
                "    lstring UnfilteredParalogy;   \"Unfiltered paralogous alignment IDs\"\n" +
                "    lstring collapsedGeneIds;   \"Collapsed Gene IDs\"\n" +
                "    lstring collapsedGeneNames;  \"Collapsed Gene Names\"\n" +
                "    string frameshift;  \"Frameshifted relative to source?\"\n" +
                "    lstring exonAnnotationSupport;   \"Exon support in reference annotation\"\n" +
                "    lstring intronAnnotationSupport;   \"Intron support in reference annotation\"\n" +
                "    string transcriptClass;    \"Transcript class\"\n" +
                "    string transcriptModes;    \"Transcript mode(s)\"\n" +
                "    string validStart;         \"Valid start codon\"\n" +
                "    string validStop;          \"Valid stop codon\"\n" +
                "    string properOrf;           \"Proper multiple of 3 ORF\"\n" +
                "    string extra_paralog;   \"Extra paralog of gene?\"\n" +
                "    )";


        BBUtils.ASTable asTable = BBUtils.parseAutosql(autosql);
        List<BBUtils.ASField> fiels = asTable.fields;
        assertEquals("name2", fiels.get(12).name);
        assertEquals("cdsStartStat", fiels.get(13).name);
        assertEquals("exonFrames", fiels.get(15).name);


    }


}
