package org.igv.tools.sort;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import org.igv.logging.*;
import org.igv.ui.IGV;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;

public class BAMSorter implements Sorter {

    private static Logger log = LogManager.getLogger(BAMSorter.class);

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
        if (this.tmpDir != null) {
            samFileWriterFactory.setTempDirectory(this.tmpDir);
        }
        if (this.maxRecords > 0) {
            samFileWriterFactory.setMaxRecordsInRam(this.maxRecords);
        }
        final SAMFileWriter writer = samFileWriterFactory.makeSAMOrBAMWriter(reader.getFileHeader(), false, outputFile);

        int count = 0;
        for (final SAMRecord rec : reader) {
            if (++count % 100000 == 0) {
                System.out.println("" + count + " records processed");   // GUI
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
        throw new RuntimeException("Not implemented");

    }

}
