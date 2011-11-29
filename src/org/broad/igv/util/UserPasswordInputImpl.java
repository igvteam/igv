package org.broad.igv.util;

import org.broad.igv.ui.IGV;
import org.broad.tribble.util.ftp.UserPasswordInput;

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
