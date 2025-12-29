package org.igv.sam.smrt;

import org.igv.sam.Alignment;
import org.igv.sam.SAMAlignment;

public class SMRTKinetics {

    private Alignment alignment;
    private short[] smrtSubreadIpd;
    private short[] smrtSubreadPw;
    private short[] smrtCcsIpd;
    private short[] smrtCcsPw;

    private short[] smrtSubreadIpdVals;
    private short[] smrtSubreadPwVals;
    private short[] smrtCcsFwdIpdVals;
    private short[] smrtCcsFwdPwVals;
    private short[] smrtCcsRevIpdVals;
    private short[] smrtCcsRevPwVals;

    public SMRTKinetics(Alignment alignment) {
        this.alignment = alignment;
    }


    /**
     * Given the compressed byte code array of CCS kinetic values from the SAM aux field do the following:
     * - Initialize frame count array, decode compressed kinetic byte codes into frame counts and copy these into it
     * - Reverse frame count sequence if reverseSequence is true
     * - If the read has a leading hard-clip, remove all frame counts from the hard-clipped region
     *
     * Returns the processed frame count array
     */
    private short[] parseSmrtKineticByteCodes(byte[] smrtKineticByteCodes, boolean reverseSequence) {
        if (smrtKineticByteCodes.length == 0) {
            return null;
        }

        int hardClipLength = alignment. getLeadingHardClipLength();
        int kineticValLength = smrtKineticByteCodes.length - hardClipLength;
        assert(kineticValLength > 0);
        short[] ccsKineticVals = new short[kineticValLength];

        int byteCodeIndex = 0;
        for (byte ccsKineticByteCode : smrtKineticByteCodes) {
            int valIndex = byteCodeIndex;
            byteCodeIndex++;
            if (reverseSequence) {
                valIndex = smrtKineticByteCodes.length - (valIndex+1);
            }
            if (valIndex < hardClipLength) continue;
            valIndex -= hardClipLength;
            ccsKineticVals[valIndex] = SMRTKineticsDecoder.lookupFrameCount(ccsKineticByteCode);
        }
        return ccsKineticVals;
    }

    /**
     * Given the original uint16 CCS kinetic frame count array from the SAM aux field do the following:
     * - Initialize the processed frame count array and copy values into it
     * - Reverse frame count sequence if reverseSequence is true
     * - If the read has a leading hard-clip, remove all frame counts from the hard-clipped region
     * <p>
     * Returns the processed frame count array
     */
    private short[] parseAuxSmrtKineticFrameCounts(short[] auxSmrtKineticFrameCounts, boolean reverseSequence) {
        if (auxSmrtKineticFrameCounts.length == 0) {
            return null;
        }

        int hardClipLength = alignment.getLeadingHardClipLength();
        int kineticValLength = auxSmrtKineticFrameCounts.length - hardClipLength;
        assert(kineticValLength > 0);
        short[] ccsKineticVals = new short[kineticValLength];

        int auxIndex = 0;
        for (short ccsKineticFrameCount : auxSmrtKineticFrameCounts) {
            int valIndex = auxIndex;
            auxIndex++;
            if (reverseSequence) {
                valIndex = auxSmrtKineticFrameCounts.length - (valIndex+1);
            }
            if (valIndex < hardClipLength) continue;
            valIndex -= hardClipLength;
            ccsKineticVals[valIndex] = ccsKineticFrameCount;
        }
        return ccsKineticVals;
    }

    private short[] getSmrtKineticsVals(short[] smrtKineticVals, String tag, boolean bytesReversedInBam) {
        if (smrtKineticVals == null && alignment.getAttribute(tag) != null) {
            final boolean reverseSequence = (bytesReversedInBam ^ alignment.isNegativeStrand());
            Object tagValue = alignment.getAttribute(tag);
            if (tagValue instanceof byte[]) {
                smrtKineticVals = parseSmrtKineticByteCodes((byte[]) (tagValue), reverseSequence);
            } else if (tagValue instanceof short[]) {
                smrtKineticVals = parseAuxSmrtKineticFrameCounts((short[]) (tagValue), reverseSequence);
            } else {
                throw new RuntimeException("Unexpected format in SMRT kinetic aux tag '" + tag + "'");
            }
        }
        return smrtKineticVals;
    }


    public short[] getSmrtSubreadIpd() {
        return getSmrtKineticsVals(smrtSubreadIpdVals, "ip", false);
    }

    public short[] getSmrtSubreadPw() {
        return getSmrtKineticsVals(smrtSubreadPwVals, "pw", false);
    }

    public short[] getSmrtCcsIpd(boolean isForwardStrand) {
        if (isForwardStrand ^ alignment.isNegativeStrand()) {
            return getSmrtKineticsVals(smrtCcsFwdIpdVals, "fi", false);
        } else {
            return getSmrtKineticsVals(smrtCcsRevIpdVals, "ri", true);
        }
    }

    public short[] getSmrtCcsPw(boolean isForwardStrand) {
        if (isForwardStrand ^ alignment.isNegativeStrand()) {
            return getSmrtKineticsVals(smrtCcsFwdPwVals, "fp", false);
        } else {
            return getSmrtKineticsVals(smrtCcsRevPwVals, "rp", true);
        }
    }
}
