/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.dev.plugin;

import javax.swing.*;

/**
 * Base class for an argument row.
 * User: jacob
 * Date: 2012-Aug-08
 */
public class ArgumentPanel extends JPanel {

    public ArgumentPanel() {
        super();
    }

    /**
     * Get the value stored for this argument. Subclasses must override
     *
     * @return
     */
    public Object getValue() {
        return null;
    }

    /**
     * Create the appropriate ArgumentPanel for this argument.
     *
     * @param argument
     * @return
     * @throws IllegalArgumentException If type not found
     */
    public static ArgumentPanel create(Argument argument) {
        ArgumentPanel panel = null;
        switch (argument.getType()) {
            case TEXT:
                panel = new TextArgument(argument);
                break;
            case FEATURE_TRACK:
                panel = new TrackArgument(argument);
                break;
            case MULTI_FEATURE_TRACK:
                panel = new MultiTrackArgument(argument);
                break;
            default:
                throw new IllegalArgumentException("Could not create ArgumentPanel for argument of type" + argument.getType());
        }
        return panel;

    }
}
