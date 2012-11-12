package org.broad.igv.hic.data;


import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.ChromosomeImpl;
import org.broad.igv.hic.HiC;
import org.broad.igv.hic.tools.Preprocessor;
import org.broad.igv.util.CompressionUtils;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.SeekableStream;

import java.io.*;
import java.util.*;

/**
 * @author jrobinso
 * @date Aug 17, 2010
 */
public class DatasetReaderV2 implements DatasetReader {

    private String path;
    private SeekableStream stream;
    private Map<String, Preprocessor.IndexEntry> masterIndex;

    private Dataset dataset = null;
    private int version = -1;
    private Map<String,int[]> fragmentSitesMap;
    private final CompressionUtils compressionUtils;

    public DatasetReaderV2(String path) throws IOException {
        this.path = path;
        this.stream = IGVSeekableStreamFactory.getStreamFor(path);
        if (this.stream != null) {
            masterIndex = new HashMap<String, Preprocessor.IndexEntry>();
            dataset = new Dataset(this);
        }
        compressionUtils = new CompressionUtils();
    }

    public static String getMagicString(String path) throws IOException {

        SeekableStream stream = null;

        try {
            stream = IGVSeekableStreamFactory.getStreamFor(path);
            LittleEndianInputStream dis = new LittleEndianInputStream(new BufferedInputStream(stream));

            String magicString = dis.readString();
            return magicString;
        } finally {
            if (stream != null) stream.close();
            ;
        }
    }

    @Override
    public String getPath() {
        return path;
    }

    @Override
    public Dataset read() throws FileNotFoundException {

        try {
            // Read the header
            LittleEndianInputStream dis = new LittleEndianInputStream(new BufferedInputStream(stream));

            String magicString = dis.readString();
            // TODO confirm magicString equals "HIC"

            version = dis.readInt();

            long masterIndexPos = dis.readLong();

            // Genome id (currently not used)
            String genomeId = dis.readString();
            dataset.setGenomeId(genomeId);


            // Read chromosome dictionary
            int nchrs = dis.readInt();
            Chromosome[] chromosomes = new Chromosome[nchrs];
            for (int i = 0; i < nchrs; i++) {
                String name = dis.readString();
                int size = dis.readInt();
                chromosomes[i] = new ChromosomeImpl(i, name, size);
            }
            dataset.setChromosomes(chromosomes);

            int nBpResolutions = dis.readInt();
            int[] bpBinSizes = new int[nBpResolutions];
            for (int i = 0; i < nBpResolutions; i++) {
                bpBinSizes[i] = dis.readInt();
            }
            dataset.setBpBinSizes(bpBinSizes);

            int nFragResolutions = dis.readInt();
            int[] fragBinSizes = new int[nFragResolutions];
            for (int i = 0; i < nFragResolutions; i++) {
                fragBinSizes[i] = dis.readInt();
            }
            dataset.setFragBinSizes(fragBinSizes);


            fragmentSitesMap = new HashMap<String, int[]>();
            if (nFragResolutions > 0) {
                for (int i = 0; i < nchrs; i++) {
                    String chr = chromosomes[i].getName();
                    int nSites = dis.readInt();
                    int [] sites = new int[nSites];
                    for (int s = 0; s < nSites; s++) {
                        sites[s] = dis.readInt();
                    }
                    fragmentSitesMap.put(chr, sites);
                }
            }

            int nHemiFragResolutions = dis.readInt();  // Reserved for future use

            // Read attribute dictionary.  Can contain arbitrary # of attributes as key-value pairs, including version

            int nAttributes = dis.readInt();
            for (int i = 0; i < nAttributes; i++) {
                String key = dis.readString();
                String value = dis.readString();
            }


            readMasterIndex(masterIndexPos, version);


        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }


        return dataset;

    }

    @Override
    public int getVersion() {
        return version;
    }

    private Map<String, Preprocessor.IndexEntry> readMasterIndex(long position, int version) throws IOException {
        stream.seek(position);
        byte[] buffer = new byte[4];
        stream.readFully(buffer);
        LittleEndianInputStream dis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
        int nBytes = dis.readInt();

        buffer = new byte[nBytes];
        stream.readFully(buffer);

        dis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
        int nEntries = dis.readInt();
        for (int i = 0; i < nEntries; i++) {
            String key = dis.readString();
            long filePosition = dis.readLong();
            int sizeInBytes = dis.readInt();
            masterIndex.put(key, new Preprocessor.IndexEntry(filePosition, sizeInBytes));
        }

        // Expected values
        Map<String, DensityFunction> expectedValuesMap = new HashMap<String, DensityFunction>();

        int nExpectedValues = dis.readInt();
        for (int i = 0; i < nExpectedValues; i++) {

            String unit = dis.readString();
            int binSize = dis.readInt();
            String key = unit + "_" + binSize;

            int nValues = dis.readInt();
            double[] values = new double[nValues];
            for (int j = 0; j < nValues; j++) {
                values[j] = dis.readDouble();
            }

            int nNormalizationFactors = dis.readInt();
            Map<Integer, Double> normFactors = new LinkedHashMap<Integer, Double>();
            for (int j = 0; j < nNormalizationFactors; j++) {
                Integer chrIdx = dis.readInt();
                Double normFactor = dis.readDouble();
                normFactors.put(chrIdx, normFactor);
            }

            DensityFunction df = new DensityFunction(unit, binSize, values, normFactors);
            expectedValuesMap.put(key, df);

        }

        dataset.setDensityFunctionMap(expectedValuesMap);

        return masterIndex;

    }


    @Override
    public Matrix readMatrix(String key) throws IOException {
        Preprocessor.IndexEntry idx = masterIndex.get(key);
        if (idx == null) {
            return null;
        }

        byte[] buffer = new byte[idx.size];
        stream.seek(idx.position);
        stream.readFully(buffer);
        LittleEndianInputStream dis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));

        int c1 = dis.readInt();
        int c2 = dis.readInt();
        Chromosome chr1 = dataset.getChromosomes()[c1];
        Chromosome chr2 = dataset.getChromosomes()[c2];

        // # of resolution levels (bp and frags)
        int nResolutions = dis.readInt();
        // dataset.setNumberZooms(nZooms);


        // TODO temporarily ignore all frag resolutions except "1f".  New UI needed for others

        List<MatrixZoomData> zdList  = new ArrayList<MatrixZoomData>();
        for (int i = 0; i < nResolutions; i++) {
            MatrixZoomData zd = new MatrixZoomData(chr1, chr2, this, dis, fragmentSitesMap);
            if(zd.getUnit() == HiC.Unit.BP || zd.getBinSize() == 1) {
                zdList.add(zd);
            }
        }
        MatrixZoomData [] zdArray =  zdList.toArray(new MatrixZoomData[zdList.size()]);

        Matrix m = new Matrix(c1, c2, zdArray);
        return m;
    }

    @Override
    public Block readBlock(int blockNumber, Preprocessor.IndexEntry idx) throws IOException {

        byte[] compressedBytes = new byte[idx.size];
        stream.seek(idx.position);
        stream.readFully(compressedBytes);

        byte[] buffer = compressionUtils.decompress(compressedBytes);
        LittleEndianInputStream dis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));

        int nRecords = dis.readInt();
        ContactRecord[] records = new ContactRecord[nRecords];
        for (int i = 0; i < nRecords; i++) {
            try {
                int bin1 = dis.readInt();
                int bin2 = dis.readInt();
                float counts = dis.readFloat();
                records[i] = new ContactRecord(blockNumber, bin1, bin2, counts);
            } catch (EOFException e) {
                nRecords = i;
                ContactRecord[] modifiedRecords = new ContactRecord[nRecords];
                System.arraycopy(records, 0, modifiedRecords, 0, nRecords);
                records = modifiedRecords;
                break;
            }

        }

        Arrays.sort(records);
        return new Block(blockNumber, records);

    }


}
