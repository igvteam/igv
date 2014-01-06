/** Copyright (c) 2007 Memorial Sloan-Kettering Cancer Center.
 **
 ** Code written by: Ethan Cerami, Benjamin Gross
 ** Authors: Ethan Cerami, Gary Bader, Chris Sander, Benjamin Gross
 **
 ** Extensively modifed by Jim Robinson for use with IGV and for supporting https.
 ** For https support you must also include "launch.html"
 **
 ** This library is free software; you can redistribute it and/or modify it
 ** under the terms of the GNU Lesser General Public License as published
 ** by the Free Software Foundation; either version 2.1 of the License, available
 ** at http://opensource.org/licenses/lgpl-2.1.php,  or
 ** any later version.
 **
 ** This library is distributed in the hope that it will be useful, but
 ** WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
 ** MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
 ** documentation provided hereunder is on an "as is" basis, and
 ** Memorial Sloan-Kettering Cancer Center
 ** has no obligations to provide maintenance, support,
 ** updates, enhancements or modifications.  In no event shall
 ** Memorial Sloan-Kettering Cancer Center
 ** be liable to any party for direct, indirect, special,
 ** incidental or consequential damages, including lost profits, arising
 ** out of the use of this software and its documentation, even if
 ** Memorial Sloan-Kettering Cancer Center
 ** has been advised of the possibility of such damage.  See
 ** the GNU Lesser General Public License for more details.
 **
 **/

/*
 * Javascript library to communication with java webstart application
 */

// globals

var igv = (function (igv) {

    igv.SCRIPT_ELEMENT_ID = "igv_script";
    var timeoutVar; // used to set/unset timeout handlers

    /**
     *
     * @param port -- the IGV port, typically 60151
     * @param command -- the command to send to the port, either "file" or "load"
     * @param paramString -- full parameter string (everything after ?),  including file to load.  See IGV user guide for
     *                       the complete list of paramters
     */
    igv.igvRequest = function (port, command, paramString) {

        var iOS = ( navigator.userAgent.match(/(iPad|iPhone|iPod)/g) ? true : false );
        if (iOS) {
            launchIGV(paramString);
            return;
        }

        //  var href = getSelfURL();
        var protocol = window.location.protocol;
        var isHTTPS = protocol == "https:";

        if (isHTTPS) {
            // We can't talk to IGV via http from this https page.  Open a new http page to do the request

            var selfURL = getSelfURL();
            var launchURL = selfURL.replace("https:", "http:").replace("webstart.js", "launchIGV.html");
            launchURL += "?port=" + port;
            launchURL += "&command=" + command;
            launchURL += "&paramString=" + encodeURIComponent(paramString);

            var salt = Math.random();
            launchURL += "&salt=" + salt;

            window.open(
                launchURL,
                "igv launch",
                "menubar=no,height=150,width=200,location=no,status=no,titlebar=no,toolbar=no",
                false);
        } else {

            var salt = Math.random(); // to prevent the browser from caching the response and preventing a relaunch if igv was shut down
            var localURL = "http://127.0.0.1:" + port + "/" + command + "?" + decodeURIComponent(paramString) + "&callback=igv.callBack();&salt=" + salt;

            //create new script
            var newScript = document.createElement("script");
            newScript.id = igv.SCRIPT_ELEMENT_ID;
            newScript.setAttribute("type", "text/javascript");
            newScript.setAttribute("src", localURL);

            // add new script to document (head section)
            var head = document.getElementsByTagName("head")[0];
            head.appendChild(newScript);
            var code = 'igv.launchIGV("' + paramString + '")';
            timeoutVar = setTimeout(code, 2000);
        }

    }

    /*
     * Handler function to launch IGV via java web start.  This handler is scheduled in the appRequest() function, and
     * is canceled by the callBack() function called in the response to the "localhost" request.  If callBack() is not
     * invoked we conclude IGV is not running and launch it via Java WebStart.
     */
    igv.launchIGV = function (queryString) {


        var protocol = window.location.protocol;
        var wsProtocol = (protocol == "https:" ? "https:" : "http:");
        var webstart_url = wsProtocol + "//www.broadinstitute.org/igv/projects/current/igv.php";

        if (queryString) {
            webstart_url += "?" + queryString;
        }

        var isChildWindow = window.location.href.indexOf("launchIGV.html") > 0;

        if (isChildWindow) {
            window.open(webstart_url);
            window.close();
        }
        else {
            window.location = webstart_url;
        }

    }

    /*
     * This function is called by IGV in the response to the GET request to load the data.  It cancels the JNLP load.
     */
    igv.callBack = function () {
        clearTimeout(timeoutVar);
        var isChildWindow = window.location.href.indexOf("launchIGV.html") > 0;
        if (isChildWindow) {
            window.close();
        }

    }


    function getSelfURL() {

        var sc = document.getElementsByTagName("script");

        for (var idx = 0; idx < sc.length; idx++) {
            var s = sc.item(idx);

            if (s.src && s.src.match(/webstart\.js$/)) {
                return s.src;
            }
        }
    }

    igv.getQueryValue = function (name, queryString) {

        if (!queryString) queryString = window.location.href;

        name = name.replace(/[\[]/, "\\\[").replace(/[\]]/, "\\\]");
        var regexS = "[\\?&]" + name + "=([^&#]*)";
        var regex = new RegExp(regexS);
        var results = regex.exec(queryString);
        if (results == null)
            return "";
        else
            return results[1];
    }

    return igv;

})(igv || {});    // Create "igv" function object if it doesn't exist


/**
 * @Deprecated.  This is kept for legacy applications, but igvRequest should be used for new development.
 *
 * This function is called from a link or button to load data into IGV.  First,  an attempt is made to load the
 * supplied data into a running IGV.  If this is not successful, as detected by a failure to cancel the timeoutHandler,
 * IGV is launched by JNLP.
 *
 * The first 2 arguments are required.  Remaining arguments are optional but must appear in the prescribed order.
 *
 * @param port -- the IGV port, typically 60151
 * @param dataUrl -- an http or ftp url to the data.
 * @param genomeID -- the genomeID,  e.g. hg18
 * @param mergeFlag -- flag to indicate if data should be merged with current IGV session, or should start a new session
 * @param locusString -- an IGV locus string, e.g. chr1:100,000-200,000  or EGFR.  See IGV doc for full details
 * @param trackName -- name for the track resulting from dataURL.  This only works for "single-track" formats, e.g. wig.
 */
function appRequest(port, dataUrl, genomeID, mergeFlag, locusString, trackName) {

    // be good and remove the previous script element
    // although, based on debugging, i'm not sure this really does anything
    var oldScript = document.getElementById(igv.SCRIPT_ELEMENT_ID);
    if (oldScript) {
        oldScript.parentNode.removeChild(oldScript);
    }

    var paramString = "file=" + dataUrl;
    if (genomeID) {
        paramString += "&genome=" + genomeID;
    }
    if (locusString) {
        paramString += "&locus=" + locusString;
    }
    if (mergeFlag) {
        paramString += "&merge=" + mergeFlag;
    }
    if (trackName) {
        paramString += "&name=" + trackName;
    }

    igv.igvRequest(port, "load", paramString);

}

