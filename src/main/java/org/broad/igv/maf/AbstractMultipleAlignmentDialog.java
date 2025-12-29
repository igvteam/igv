package org.broad.igv.maf;

import java.awt.*;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 2/15/13
 *         Time: 7:24 PM
 */
public abstract class AbstractMultipleAlignmentDialog extends org.broad.igv.ui.IGVDialog  {

    public AbstractMultipleAlignmentDialog(Frame parent, boolean modal) {
        super(parent, modal);
    }

    abstract List<String> getSelectedSpecies();

    public abstract boolean isCancelled();
}
