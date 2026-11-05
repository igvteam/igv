package org.igv.ucsc.bb.codecs;

import org.igv.feature.BasicFeature;
import org.igv.feature.IGVFeature;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.ucsc.bb.BBUtils;
import org.igv.ucsc.bb.BedData;

import java.util.List;

/**
 * Base class for codecs that decode a bigBed record directly into a specialized feature type.  Provides the autosql
 * field lookup helpers and the fallback for malformed records: a plain bed feature, with a warning logged for the
 * first few.
 */
abstract class BBRecordCodec implements BBCodec {

    private static Logger log = LogManager.getLogger(BBRecordCodec.class);

    private static final int MAX_WARNINGS = 10;

    private final String label;
    private final BBBedCodec bedCodec;
    private int warningCount = 0;

    BBRecordCodec(String label, int standardFieldCount, BBUtils.ASTable autosql) {
        this.label = label;
        this.bedCodec = new BBBedCodec(standardFieldCount, autosql);
    }

    /**
     * Decode the record, or return null if it cannot be decoded.  NumberFormatExceptions are treated as null.
     */
    abstract IGVFeature decodeRecord(BedData data);

    @Override
    public IGVFeature decode(BedData data) {
        try {
            IGVFeature record = decodeRecord(data);
            if (record != null) {
                return record;
            }
        } catch (NumberFormatException e) {
            // Fall through
        }
        if (warningCount++ < MAX_WARNINGS) {
            log.warn("Malformed " + label + " record at " + data.chr + ":" + data.chromStart + "-" + data.chromEnd);
        }
        try {
            return bedCodec.decode(data);
        } catch (Exception e) {
            return new BasicFeature(data.chr, data.chromStart, data.chromEnd);
        }
    }

    /**
     * Return the index of the named autosql field in "restOfFields" (which starts at the 4th bed column), or -1 if
     * not defined.
     */
    static int fieldIndex(BBUtils.ASTable autosql, String name) {
        List<BBUtils.ASField> fields = autosql.fields;
        for (int i = 3; i < fields.size(); i++) {
            if (name.equals(fields.get(i).name)) {
                return i - 3;
            }
        }
        return -1;
    }

    static String stringValue(String[] tokens, int idx, String defaultValue) {
        return idx < 0 || idx >= tokens.length || tokens[idx].isEmpty() ? defaultValue : tokens[idx];
    }

    static int intValue(String[] tokens, int idx) {
        return Integer.parseInt(stringValue(tokens, idx, "0"));
    }
}
