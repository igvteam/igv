package org.broad.igv.sam.cram;

import com.google.common.io.LittleEndianDataInputStream;
import com.google.common.io.LittleEndianDataOutputStream;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.util.HttpUtils;

import java.io.*;
import java.util.Arrays;
import java.util.Comparator;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.zip.*;

/**
 * Created by jrobinso on 6/25/17.
 * <p>
 * Some static methods for managing the CRAM reference sequence cache
 */

public class ReferenceDiskCache {

    private static Logger log = Logger.getLogger(ReferenceDiskCache.class);

    private static final ExecutorService threadExecutor = Executors.newFixedThreadPool(1);

    public static void saveSequence(final String genomeId, final String chr, final byte[] bytes) throws IOException {

        threadExecutor.submit(() -> {
            {
                File cacheDir = getCacheDirectory();

                FileOutputStream fos = null;
                final File outputFile = new File(cacheDir, getFileName(genomeId, chr));
                try {

                    final FileOutputStream out = new FileOutputStream(outputFile);

                    LittleEndianDataOutputStream dos = new LittleEndianDataOutputStream(out);
                    dos.writeInt(bytes.length);

                    DeflaterOutputStream gzipOutputStream = new DeflaterOutputStream(out);
                    gzipOutputStream.write(bytes);
                    gzipOutputStream.flush();
                    gzipOutputStream.close();

                    checkCacheSize();

                } catch (Exception e) {

                    log.error("Error saving CRAM reference sequence", e);
                    outputFile.delete();

                } finally {
                    if (fos != null) {
                        try {
                            fos.close();
                        } catch (IOException e) {

                        }
                    }
                }
            }
        });
    }

    public static byte[] readSequence(String genomeId, String chr) throws IOException {

        File cacheDirectory = getCacheDirectory();
        File seqFile = new File(cacheDirectory, getFileName(genomeId, chr));

        if (!seqFile.exists()) return null;

        FileInputStream fis = null;

        try {
            fis = new FileInputStream(seqFile);
            LittleEndianDataInputStream dis = new LittleEndianDataInputStream(fis);

            int size = dis.readInt();

            InflaterInputStream gis = new InflaterInputStream(fis);   //new ByteArrayInputStream(compressedBytes));
            byte[] buffer = new byte[size];
            HttpUtils.readFully(gis, buffer);
            return buffer;
        } finally {
            if (fis != null) fis.close();
        }


    }

    private static File getCacheDirectory() {

        String rootDirectoryString = PreferencesManager.getPreferences().get(Constants.CRAM_CACHE_DIRECTORY);

        File rootDirectory;
        if (rootDirectoryString != null) {
            rootDirectory = new File(rootDirectoryString);
        } else {
            rootDirectory = new File(DirectoryManager.getIgvDirectory(), "cram");
        }

        if (!rootDirectory.exists()) {
            rootDirectory.mkdir();
        }

        return rootDirectory;
    }


    private static synchronized void checkCacheSize() {

        File cacheDir = getCacheDirectory();
        if(!cacheDir.exists()) {
            return;
        }

        File[] files = cacheDir.listFiles((dir, name) -> {return name.toLowerCase().endsWith(".bin");});
        Arrays.sort(files, Comparator.comparingLong(File::lastModified).reversed());

        int maxSize = PreferencesManager.getPreferences().getAsInt(Constants.CRAM_CACHE_SIZE) * 1000;
        int totalSize = 0;
        for(File f : files) {
            if(totalSize > maxSize) {
                f.delete();
            }
            else {
                totalSize += f.length();
            }
        }
    }


    public static void deleteCache(String genomeId, String chr) {
        (new File(getFileName(genomeId, chr))).delete();
    }


    private static String getFileName(String genomeId, String chr) {
        // genomeIds can be full paths and other illegal filename strings.
        return String.valueOf(genomeId.hashCode()) + "-" + chr + ".bin";
    }
}
