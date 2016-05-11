/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 University of California San Diego
 * Author: Jim Robinson
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

package org.broad.igv.tools.sort;

import com.mysql.jdbc.NotImplemented;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import org.apache.log4j.Logger;
import org.broad.igv.ui.IGV;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;

public class BAMSorter implements Sorter {

    private static Logger log = Logger.getLogger(BAMSorter.class);

    File inputFile;
    File outputFile;
    private File tmpDir;
    private int maxRecords = -1;

    public BAMSorter(File inputFile, File outputFile) {
        this.inputFile = inputFile;
        this.outputFile = outputFile;
    }

    @Override
    public void run() throws IOException {

        SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
        final SamReader reader = SamReaderFactory.makeDefault().open(inputFile);


        reader.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.coordinate);
        SAMFileWriterFactory samFileWriterFactory = new SAMFileWriterFactory();
        if(this.tmpDir != null) {
            samFileWriterFactory.setTempDirectory(this.tmpDir);
        }
        if(this.maxRecords > 0) {
            samFileWriterFactory.setMaxRecordsInRam(this.maxRecords);
        }
        final SAMFileWriter writer = samFileWriterFactory.makeSAMOrBAMWriter(reader.getFileHeader(), false, outputFile);

        int count = 0;
        for (final SAMRecord rec : reader) {
            if(++count % 100000 == 0) {
                if(IGV.hasInstance()) {
                    System.out.println("" + count + " records processed");   // GUI
                }
                else {
                    log.info("" + count + " records processed");   // Command line
                }
            }
            writer.addAlignment(rec);
        }

        CloserUtil.close(reader);
        writer.close();
    }


    @Override
    public void setTmpDir(File tmpDir) {

        this.tmpDir = tmpDir;

    }

    @Override
    public void setMaxRecords(int maxRecords) {
        this.maxRecords = maxRecords;
    }

    @Override
    public void setComparator(Comparator<SortableRecord> comparator) {
        throw new NotImplementedException();

    }

}
