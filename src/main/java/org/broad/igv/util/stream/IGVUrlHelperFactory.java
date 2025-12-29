package org.broad.igv.util.stream;

import htsjdk.tribble.util.URLHelper;
//import htsjdk.tribble.util.URLHelperFactory;
import htsjdk.tribble.util.URLHelperFactory;
import org.broad.igv.logging.*;
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

