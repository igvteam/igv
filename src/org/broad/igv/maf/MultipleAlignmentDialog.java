package org.broad.igv.maf;

import java.awt.*;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 2/15/13
 *         Time: 8:39 PM
 */
public class MultipleAlignmentDialog extends AbstractMultipleAlignmentDialog {


    public MultipleAlignmentDialog(Frame parent, boolean modal) {
        super(parent, modal);
    }

    @Override
    List<String> getSelectedSpecies() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean isCancelled() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
