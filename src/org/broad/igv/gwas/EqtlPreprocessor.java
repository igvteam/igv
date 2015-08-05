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

package org.broad.igv.gwas;

import org.broad.igv.tdf.BufferedByteWriter;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.tools.sort.Sorter;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.CompressionUtils;

import java.io.*;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

/**
 * Class for experimenting with binary EQTL formats.
 * <p/>
 * SNP	SNP_Chr	SNP_Pos	Gen_ID	Gene_Name	Gene_Pos	T_Stat	P_Val	Q_Val
 * rs3094317	1	739528	ENSG00000237491.2	RP11-206L10.6	714162	-4.41310870534187	3.457968769964e-05	0.0254439762182497
 *
 SNP	SNP_Chr	SNP_Pos	Gen_ID	Gene_Name	Gene_Pos	T_Stat	Beta	P_Val	min(p)	EmpP	nom_thresh
 chr2:202672143:I	2	202672143	ENSG00000003393.10	ALS2	202645912	-4.90641599377695	-0.410144871510459	4.16077158322729e-06	4.1199296E-8	9 *
 */
public class EqtlPreprocessor {

    static int version = 0;

    long indexPositionPosition;
    int bytesWritten;
    BufferedByteWriter currentChrBuffer;
    String currentChr;
    private FileOutputStream fos;
    private Map<String, IndexEntry> chrPositionMap;
    private CompressionUtils compressionUtils;


    static class IndexEntry {

        IndexEntry(long position, int size) {
            this.position = position;
            this.size = size;
        }

        long position;
        int size;
    }

    public static void main(String[] args) throws IOException {

        File[] files = (new File(args[0])).listFiles();
        for (File file : files) {
            if (file.getName().endsWith(".eqtl")) {

                File sortedFile = new File(file.getAbsolutePath() + ".sorted.eqtl");
                Sorter sorter = Sorter.getSorter(file, sortedFile);
                sorter.run();

                (new EqtlPreprocessor()).process(sortedFile.getAbsolutePath(), file.getAbsolutePath() + ".bin");
            }
        }
    }

    public void process(String inputFile, String outputFile) throws IOException {

        bytesWritten = 0;
        BufferedReader br = null;
        fos = null;
        currentChr = null;
        currentChrBuffer = null;
        chrPositionMap = new HashMap<String, IndexEntry>();
        compressionUtils = new CompressionUtils();
        EQTLCodec codec = new EQTLCodec(null);

        try {
            fos = new FileOutputStream(outputFile);
            writeHeader();

            br = new BufferedReader(new FileReader(inputFile));
            String nextLine = br.readLine();

            while ((nextLine = br.readLine()) != null) {

                EQTLFeature feature = codec.decode(nextLine);
                String chr = feature.getChr();
                if (!chr.equals(currentChr)) {
                    if (currentChrBuffer != null) {
                        System.out.println(currentChr);
                        writeChrBuffer();
                    }
                    currentChr = chr;
                    currentChrBuffer = new BufferedByteWriter();
                }

                currentChrBuffer.put(feature.encodeBinary());

            }
            if (currentChrBuffer != null) {
                writeChrBuffer();
            }

            long position = bytesWritten;
            int nBytes = writeIndex();
            fos.close();

            writeIndexPosition(outputFile, position, nBytes);

        } finally {
            if (br != null) br.close();
        }

    }

    private void writeChrBuffer() throws IOException {
        // compress and write
        byte[] rawBytes = currentChrBuffer.getBytes();
        byte[] compressedBytes = rawBytes; //compressionUtils.compress(rawBytes);
        int position = bytesWritten;
        int size = compressedBytes.length;
        write(compressedBytes);
        chrPositionMap.put(currentChr, new IndexEntry(position, size));

    }

    private int writeIndex() throws IOException {
        // compress and write

        BufferedByteWriter buff = new BufferedByteWriter();
        buff.putInt(chrPositionMap.size());
        for (Map.Entry<String, IndexEntry> entry : chrPositionMap.entrySet()) {
            IndexEntry ie = entry.getValue();
            buff.putNullTerminatedString(entry.getKey());
            buff.putLong(ie.position);
            buff.putInt(ie.size);
        }

        byte[] bytes = buff.getBytes();
        write(bytes);
        return bytes.length;

    }

    private void writeHeader() throws IOException {

        // Magic number -- 4 bytes
        byte[] magicNumber = new byte[]{'E', 'Q', 'T', 'L'};

        BufferedByteWriter buffer = new BufferedByteWriter();
        buffer.put(magicNumber);
        buffer.putInt(version);
        // Reserve space for the master index pointer and byte count.
        // The actual values will be written at the end
        indexPositionPosition = buffer.bytesWritten();
        buffer.putLong(0l);   // File position for start of index
        buffer.putInt(0);     // Size in bytes of index

        final byte[] bytes = buffer.getBytes();
        write(bytes);
    }

    private void writeIndexPosition(String file, long indexPosition, int nbytes) {
        try {
            RandomAccessFile raf = new RandomAccessFile(file, "rw");
            raf.getChannel().position(indexPositionPosition);

            // Write as little endian
            BufferedByteWriter buffer = new BufferedByteWriter();
            buffer.putLong(indexPosition);
            buffer.putInt(nbytes);
            raf.write(buffer.getBytes());
            raf.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }


    private void write(byte[] bytes) throws IOException {
        fos.write(bytes);
        bytesWritten += bytes.length;
    }

}
