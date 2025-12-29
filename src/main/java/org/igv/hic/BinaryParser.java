package org.igv.hic;

    import java.nio.ByteBuffer;
    import java.nio.ByteOrder;
    import java.nio.charset.StandardCharsets;

    /**
     * A parser for reading binary data from a ByteBuffer, converted from a JavaScript utility.
     * It assumes little-endian byte order by default, which can be overridden in the constructor.
     */
    public class BinaryParser {

        private final ByteBuffer buffer;

        public BinaryParser(ByteBuffer buffer) {
            this(buffer, true);
        }

        public BinaryParser(ByteBuffer buffer, boolean littleEndian) {
            this.buffer = buffer;
            if (littleEndian) {
                this.buffer.order(ByteOrder.LITTLE_ENDIAN);
            } else {
                this.buffer.order(ByteOrder.BIG_ENDIAN);
            }
        }

        public int available() {
            return buffer.remaining();
        }

        public boolean hasNext() {
            return buffer.hasRemaining();
        }

        /**
         * Reads an unsigned byte and returns it as an int.
         */
        public int getByte() {
            return buffer.get() & 0xFF;
        }

        public short getShort() {
            return buffer.getShort();
        }

        /**
         * Reads an unsigned short and returns it as an int.
         */
        public int getUShort() {
            return buffer.getShort() & 0xFFFF;
        }

        public int getInt() {
            return buffer.getInt();
        }

        /**
         * Reads an unsigned int and returns it as a long.
         */
        public long getUInt() {
            return buffer.getInt() & 0xFFFFFFFFL;
        }

        public long getLong() {
            return buffer.getLong();
        }

        public float getFloat() {
            return buffer.getFloat();
        }

        public double getDouble() {
            return buffer.getDouble();
        }

        /**
         * Reads a null-terminated string.
         */
        public String getString() {
            return getString(-1);
        }

        /**
         * Reads a string, either null-terminated or of a fixed length.
         *
         * @param len The fixed length. If -1, reads until a null terminator.
         */
        public String getString(int len) {
            StringBuilder sb = new StringBuilder();
            int i = 0;
            while (buffer.hasRemaining() && (len < 0 || i < len)) {
                byte c = buffer.get();
                if (c == 0) {
                    break;
                }
                sb.append((char) c);
                if (len > 0) {
                    i++;
                }
            }
            return sb.toString();
        }

        public String getFixedLengthString(int len) {
            byte[] bytes = new byte[len];
            buffer.get(bytes);
            int firstNull = -1;
            for (int i = 0; i < bytes.length; i++) {
                if (bytes[i] == 0) {
                    firstNull = i;
                    break;
                }
            }
            if (firstNull != -1) {
                return new String(bytes, 0, firstNull, StandardCharsets.UTF_8);
            } else {
                return new String(bytes, StandardCharsets.UTF_8);
            }
        }

        public String getFixedLengthTrimmedString(int len) {
            byte[] bytes = new byte[len];
            buffer.get(bytes);
            return new String(bytes, StandardCharsets.UTF_8).trim();
        }

        public void skip(int n) {
            buffer.position(buffer.position() + n);
        }

        /**
         * Return a bgzip (bam and tabix) virtual pointer.
         */
        public VPointer getVPointer() {
            long rawValue = buffer.getLong();
            long block = rawValue >>> 16;
            int offset = (int) (rawValue & 0xFFFF);
            return new VPointer(block, offset);
        }

        public int position() {
            return buffer.position();
        }

        /**
         * Represents a BGZF virtual pointer with a block address and an offset within the block.
         */
        public static class VPointer {
            public final long block;
            public final int offset;

            public VPointer(long block, int offset) {
                this.block = block;
                this.offset = offset;
            }

            public boolean isLessThan(VPointer vp) {
                return this.block < vp.block || (this.block == vp.block && this.offset < vp.offset);
            }

            public boolean isGreaterThan(VPointer vp) {
                return this.block > vp.block || (this.block == vp.block && this.offset > vp.offset);
            }

            @Override
            public String toString() {
                return "" + this.block + ":" + this.offset;
            }
        }
    }