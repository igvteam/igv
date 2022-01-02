package org.broad.igv.ext.annotate;

import org.broad.igv.ext.IExtension;
import org.broad.igv.sam.AlignmentBlock;
import org.broad.igv.sam.SAMAlignment;

public interface IAlignmentBlockAnnotationExtension extends IExtension {

    public void appendBlockQualityAnnotation(SAMAlignment samAlignment, AlignmentBlock block, StringBuffer buf);
    public void appendBlockAttrAnnotation(SAMAlignment samAlignment, AlignmentBlock block, int offset, StringBuffer buf);
}