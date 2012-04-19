/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package edu.mit.broad.prodinfo.multiplealignment;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Stack;

import edu.mit.broad.prodinfo.datastrutures.IntervalTree;
import edu.mit.broad.prodinfo.genomicplot.ParseException;
import edu.mit.broad.prodinfo.multiplealignment.MAFAlignment.MAFHeader;
import edu.mit.broad.prodinfo.multiplealignment.MAFAlignment.MAFMultipleAlignmentBlock;
import edu.mit.broad.prodinfo.multiplealignment.MultipleAlignment.AlignedSequence;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.SeekableStreamFactory;

public class MAFIO implements MultipleAlignmentIO {

    SeekableStream fileHandle;
    String alignmentFile;
    private IntervalTree<Long> index;

    public MAFIO(String alignmentFile, boolean keepFileHandle) throws IOException, ParseException {
        // alignment = createUnloadedAlignment(alignmentFile);
        this.alignmentFile = alignmentFile;
        if (keepFileHandle) {
            fileHandle = SeekableStreamFactory.getStreamFor(alignmentFile);
        }
        checkIndex(alignmentFile);
    }

    public MAFIO() {
        super();
    }

    public String getPreferredFileExtension() {
        return ".maf";
    }

    public void destroyFileHandle() throws IOException {
        if (fileHandle != null) {
            fileHandle.close();
        }
    }

    public MAFAlignment load(List<String> sequencesToLoad, int start, int end) throws IOException, ParseException {
        //RandomAccessFile handle = fileHandle;
        MAFAlignment alignment = null;
        boolean closeFile = false;
        if (fileHandle == null) { // It would be nice to test whether the handle is closed if not null....
            closeFile = true;
            fileHandle = SeekableStreamFactory.getStreamFor(alignmentFile);
        }
        try {
            alignment = new MAFAlignment(index);
            alignment.load(fileHandle, start, end, sequencesToLoad);
        } finally {
            if (closeFile) {
                fileHandle.close();
                fileHandle = null;
            }
        }
        return alignment;
    }

    public void write(BufferedWriter bw, MultipleAlignment ma)
            throws IOException {
        // TODO Auto-generated method stub

    }

    public void write(BufferedWriter bw, MultipleAlignment ma, List<String> sequenceOrder)
            throws IOException {
        // TODO Auto-generated method stub
    }


    private void checkIndex(String fileName) throws IOException, ParseException {
        String idxFileName = fileName + ".index";
        if (FileUtils.resourceExists(idxFileName)) {
            loadIndex(idxFileName);
        } else {
            if (FileUtils.isRemote(fileName)) {
                throw new RuntimeException("Index file not found: " + idxFileName);
            } else {
                createIndex(fileName);
                writeIndex(idxFileName);
            }
        }
    }

    /**
     * Load every 100 lines
     *
     * @param idxFile
     * @throws IOException
     */
    public void loadIndex(String idxFile) throws IOException {
        index = new IntervalTree<Long>();
        BufferedReader br = null;
        try {
            br = ParsingUtils.openBufferedReader(idxFile);
            //System.err.println("Loading index " + idxFile);
            String line = null;
            int l = 0;
            while ((line = br.readLine()) != null) {
                if (l % 100 == 0) {
                    String[] info = line.split("\t");
                    int start = Integer.parseInt(info[0]);
                    int end = Integer.parseInt(info[1]) + start;
                    long offset = Long.parseLong(info[2]);
                    //System.out.println("Read line "+ l++);
                    index.put(start, end, offset);
                }

                l++;
            }
        } finally {
            if (br != null) br.close();
        }

    }

    /**
     * Create an index for the MAF file.
     * <p/>
     * Example MAF lines:
     * a score=34237.000000
     * s hg19.chr1     10917 479 + 249250621 gagaggc
     * s panTro2.chr15 13606 455 - 100063
     *
     * @param alignmentFile
     * @throws IOException
     */
    public void createIndex(String alignmentFile) throws IOException {
        RandomAccessFile raf = new RandomAccessFile(alignmentFile, "r");
        index = new IntervalTree<Long>();//LinkedHashMap<Integer, Long>();
        String line;
        String[] lineInfo = null;
        long lastOffset = 0;
        try {
            boolean readNext = false;
            while ((line = raf.readLine()) != null) {
                //Ignore all comment lines
                if (line.startsWith("#") || line.trim().length() == 0) {
                    continue;
                }

                if (line.startsWith("a ")) {
                    readNext = true;
                    lastOffset = raf.getFilePointer() - line.getBytes().length - 1;
                } else if (line.startsWith("s ")) {
                    if (readNext) {
                        lineInfo = line.split("\\s+");
                        int start = Integer.parseInt(lineInfo[2]);
                        int end = Integer.parseInt(lineInfo[3]) + start;
                        index.put(start, end, lastOffset);
                    }
                    readNext = false;
                } else if (line.startsWith("i ")) {
                    //We do not handle information lines yet.
                    continue;
                } else if (line.startsWith("q ")) {
                    //We do not handle quality lines yet.
                    continue;
                } else {
                    readNext = false;
                }

            }
        } finally {
            try {
                raf.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

    }


    public void writeIndex(String indexFileName) throws IOException {
        Iterator<IntervalTree.Node<Long>> idxEntryIt = index.iterator();
        BufferedWriter bw = new BufferedWriter(new FileWriter(indexFileName));
        while (idxEntryIt.hasNext()) {
            IntervalTree.Node<Long> entry = idxEntryIt.next();
            bw.write(String.valueOf(entry.getStart()));
            bw.write("\t" + String.valueOf(entry.getEnd() - entry.getStart()));
            bw.write("\t" + String.valueOf(entry.getValue()));
            bw.newLine();
        }
        bw.close();
    }


}
