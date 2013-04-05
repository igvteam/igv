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

package org.broad.igv.util;

import net.sf.samtools.seekablestream.UserPasswordInput;
import org.broad.igv.ui.IGV;

/**
 * @author Jim Robinson
 * @date 10/4/11
 */
public class UserPasswordInputImpl implements UserPasswordInput {

    String host;
    String password;
    String user;

    public void setHost(String host) {
        this.host = host;
    }

    public boolean showDialog() {

        UserPasswordDialog dlg = new UserPasswordDialog(IGV.getMainFrame(), user, host);
        dlg.setVisible(true);

        if (dlg.isCanceled()) {
            dlg.dispose();
            return false;
        } else {
            user = dlg.getUser();
            password = dlg.getPassword();
            dlg.dispose();
            return true;
        }

    }

    public String getUser() {
        return user;
    }

    public String getPassword() {
        return password;
    }
}
