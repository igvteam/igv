package org.broad.igv.sam;

/// Logic to decode SMRT kinetic values from uint8 compressed code to an approximation of the
/// original frame count.
///
public class SMRTKineticsDecoder {

    private static final int base = 1 << 6;

    private final short[] frameCounts;

    public SMRTKineticsDecoder() {
        frameCounts = new short[256];
        for (int i = 0; i < 256; ++i) {
            frameCounts[i] = decodeFrameCount((byte) i);
        }
    }

    /// Convert compressed byte to frame count from cached LUT
    public short lookupFrameCount(byte val) {
        return frameCounts[val & 0xFF];
    }

    /// Compute frame count from compressed byte
    ///
    /// 8-bit compressed V1 codec:
    ///
    ///   xxmm'mmmm
    ///
    /// where
    /// * 'x' represents the exponent (2 MSBs)
    /// * 'm' represents the mantissa (6 LSBs)
    ///
    private static short decodeFrameCount(byte val) {
        int mantissa = val & 0b0011_1111;
        int exponent = (val & 0b1100_0000) >> 6;
        int exponentVal = 1 << exponent;
        return (short) (base * (exponentVal - 1) + exponentVal * mantissa);
    }
}
