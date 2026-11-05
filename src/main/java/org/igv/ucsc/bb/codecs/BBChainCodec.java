package org.igv.ucsc.bb.codecs;

import org.igv.Globals;
import org.igv.feature.PSLRecord;
import org.igv.feature.Strand;
import org.igv.ucsc.bb.BBUtils;
import org.igv.ucsc.bb.BedData;

/**
 * Codec for UCSC "bigChain" files (bed6+6).  Records are decoded into PSLRecords, with the query sequence name as
 * the feature name and the chain id as the identifier.  A bigChain record has no blocks, these are supplied by the
 * companion "bigLink" file (see BBChainFeatureSource).
 */
public class BBChainCodec extends BBRecordCodec {

    // Indices into "restOfFields" (which starts at the 4th bed column) for the standard bed6 fields
    private static final int NAME = 0, SCORE = 1, STRAND = 2;

    private final int tSize, qName, qSize, qStart, qEnd, chainScore;

    public BBChainCodec(int standardFieldCount, BBUtils.ASTable autosql) {
        super("bigChain", standardFieldCount, autosql);
        tSize = fieldIndex(autosql, "tSize");
        qName = fieldIndex(autosql, "qName");
        qSize = fieldIndex(autosql, "qSize");
        qStart = fieldIndex(autosql, "qStart");
        qEnd = fieldIndex(autosql, "qEnd");
        chainScore = fieldIndex(autosql, "chainScore");
    }

    @Override
    PSLRecord decodeRecord(BedData data) {

        String[] tokens = Globals.tabPattern.split(data.restOfFields, -1);
        String queryName = stringValue(tokens, qName, null);
        if (tokens.length <= STRAND || queryName == null) {
            return null;
        }

        String chainId = tokens[NAME];
        int queryStart = intValue(tokens, qStart);
        int queryEnd = intValue(tokens, qEnd);
        int querySize = intValue(tokens, qSize);
        Strand strand = tokens[STRAND].startsWith("+") ? Strand.POSITIVE : Strand.NEGATIVE;

        PSLRecord f = new PSLRecord();
        f.setChr(data.chr);
        f.setStart(data.chromStart);
        f.setEnd(data.chromEnd);
        f.setName(queryName);
        f.setIdentifier(chainId);
        f.setStrand(strand);
        f.setScore(Float.parseFloat(tokens[SCORE]));
        f.setTSize(intValue(tokens, tSize));
        f.setQSize(querySize);
        f.setqStart(queryStart);
        f.setqEnd(queryEnd);

        StringBuffer desc = new StringBuffer();
        desc.append("chain id = " + chainId);
        desc.append("<br>");
        desc.append("chain score = " + stringValue(tokens, chainScore, ""));
        desc.append("<br>");
        // Query coordinates of '-' strand chains are on the reverse strand, convert for display as UCSC does
        if (strand == Strand.NEGATIVE && querySize > 0) {
            desc.append("query = " + queryName + ":" + (querySize - queryEnd + 1) + "-" + (querySize - queryStart) + " (-)");
        } else {
            desc.append("query = " + queryName + ":" + (queryStart + 1) + "-" + queryEnd);
        }
        f.setDescription(desc.toString());

        return f;
    }
}
