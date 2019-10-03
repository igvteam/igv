/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.util.stream;

import htsjdk.tribble.util.URLHelper;
//import htsjdk.tribble.util.URLHelperFactory;
import htsjdk.tribble.util.URLHelperFactory;
import org.apache.log4j.Logger;
import org.broad.igv.util.HttpUtils;

import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 * @date Jul 6, 2011
 */
public class IGVUrlHelperFactory implements URLHelperFactory {

    private static IGVUrlHelperFactory theInstance;

    public static IGVUrlHelperFactory getInstance() {
        if(theInstance == null) {
            theInstance = new IGVUrlHelperFactory();
        }
        return theInstance;
    }

    private IGVUrlHelperFactory() {}

    //@Override
    public URLHelper getHelper(URL url) {
        return new IGVUrlHelper(url);
    }


}

