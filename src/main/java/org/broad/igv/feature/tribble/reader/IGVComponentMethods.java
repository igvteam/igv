package org.broad.igv.feature.tribble.reader;

import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.util.ParsingUtils;

import java.io.IOException;

public class IGVComponentMethods extends AbstractFeatureReader.ComponentMethods {


    private static  IGVComponentMethods instance;

    public static IGVComponentMethods getInstance() {
        if(instance == null) {
            instance = new IGVComponentMethods();
        }
        return instance;
    }

    private IGVComponentMethods() {}

    public boolean isTabix(String resourcePath, String indexPath) throws IOException {

        boolean isRemote = (resourcePath.startsWith("gs://") || resourcePath.startsWith("http://") || resourcePath.startsWith("https://"));

        if(isRemote) {
            return indexPath != null && IOUtil.hasBlockCompressedExtension(resourcePath);
        }
        else {
            if (indexPath == null & !isRemote) {
                indexPath = ParsingUtils.appendToPath(resourcePath, FileExtensions.TABIX_INDEX);
            }
            return indexPath != null && IOUtil.hasBlockCompressedExtension(resourcePath) && ParsingUtils.resourceExists(indexPath);
        }
    }

}
