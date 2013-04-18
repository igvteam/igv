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

package org.broad.igv.util.stream;

import org.broad.igv.util.HttpUtils;
import org.broad.tribble.util.URLHelper;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 * @date Jul 6, 2011
 */
public class IGVUrlHelper implements URLHelper {

    URL url;

    public IGVUrlHelper(URL url) {
        this.url = url;
    }

    public URL getUrl() {
        return url;
    }

    public long getContentLength() throws IOException {
        return HttpUtils.getInstance().getContentLength(url);
    }

    public InputStream openInputStream() throws IOException {
        return HttpUtils.getInstance().openConnectionStream(url);
    }

    public InputStream openInputStreamForRange(long start, long end) throws IOException {

        String byteRange = "bytes=" + start + "-" + end;
        Map<String, String> params = new HashMap();
        params.put("Range", byteRange);
        return HttpUtils.getInstance().openConnectionStream(url, params);
    }

    public boolean exists() {
        return HttpUtils.getInstance().resourceAvailable(url);
    }
}
