/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.gs;

import org.junit.Ignore;

import java.net.Authenticator;
import java.net.PasswordAuthentication;

/**
 * @author jacob
 * @date 2013-May-15
 */
@Ignore
public class GSTestAuthenticator extends Authenticator {

    @Override
    protected PasswordAuthentication getPasswordAuthentication() {
        return new PasswordAuthentication("igvtest", "igvtest".toCharArray());
    }

    public static void setTokenAuthentication(){
        Authenticator.setDefault(null);
        GSUtils.setGSUser("igvtest");
        GSUtils.setGSToken("HnR9rBShNO4dTXk8cKXVJT98Oe0jWVY+");
    }
}
