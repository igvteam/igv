// $Id: webstart.js,v 1.13 2008-01-14 15:05:28 grossben Exp $
//------------------------------------------------------------------------------
/** Copyright (c) 2007 Memorial Sloan-Kettering Cancer Center.
 **
 ** Code written by: Ethan Cerami, Benjamin Gross
 ** Authors: Ethan Cerami, Gary Bader, Chris Sander, Benjamin Gross
 ** Modified by Jim Robinson for use with IGV
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
var SCRIPT_ELEMENT_ID = "igv";
var timeoutVar; // used to set/unset timeout handlers
var aWindow;
var origProtocol;


/*
 * This function is called by IGV in the response to the GET request to load the data.  It cancels the JNLP load.
 */
function callBack() {
    clearTimeout(timeoutVar);
}

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

    // be good and remove the previous cytoscape script element
    // although, based on debugging, i'm not sure this really does anything
    var oldScript = document.getElementById(SCRIPT_ELEMENT_ID);
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

    igvRequest(port, "load", paramString);

}

/**
 *
 * @param port -- the IGV port, typically 60151
 * @param command -- the command to send to the port, either "file" or "load"
 * @param paramString -- full parameter string (everything after ?),  including file to load.  See IGV user guide for
 *                       the complete list of paramters
 */
function igvRequest(port, command, paramString, origProtocol) {

    var protocol = window.location.protocol;
    if(!origProtocol) origProtocol = protocol;

    if (aWindow) {
        aWindow.close();
    }

    var isIpad = navigator.userAgent.match(/iPad/i) != null;
    if (isIpad) {
        launchIGV(paramString);
        return;
    }

    //  var href = getSelfURL();
    var isHTTPS = protocol == "https:";

    if (isHTTPS) {
        // We can't talk to IGV via http from this https page.  Open a new http page to do the request

        var selfURL = getSelfURL();
        var launchURL = selfURL.replace("https:", "http:").replace("webstart.js", "launch.html");
        launchURL += "?port=" + port;
        launchURL += "&command=" + command;
        launchURL += "&protocol=" + protocol;
        launchURL += "&paramString=" + encodeURIComponent(paramString);

        var salt = Math.random();
        launchURL += "&salt=" + salt;

        aWindow = window.open(
            launchURL,
            "igv launch",
            "menubar=no,height=150,width=200,location=no,status=no,titlebar=no,toolbar=no",
            false);
//        aWindow.document.close();
    }
    else {

        var salt = Math.random(); // to prevent the browser from caching the response and preventing a relaunch if igv was shut down
        var localURL = "http://127.0.0.1:" + port + "/" + command + "?paramString=" + decodeURIComponent(paramString) + "&callback=callBack();&salt=" + salt;

        //create new script
        var newScript = document.createElement("script");
        newScript.id = SCRIPT_ELEMENT_ID;
        newScript.setAttribute("type", "text/javascript");
        newScript.setAttribute("src", localURL);

        // add new script to document (head section)
        var head = document.getElementsByTagName("head")[0];
        head.appendChild(newScript);
        timeoutVar = setTimeout('launchIGV("' + origProtocol + '", "' + paramString + '")', 2000);
    }

}

/*
 * Handler function to launch IGV via java web start.  This handler is scheduled in the appRequest() function, and
 * is canceled by the callBack() function called in the response to the "localhost" request.  If callBack() is not
 * invoked we conclude IGV is not running and launch it via Java WebStart.
 */
function launchIGV(protocol, queryString) {


    var webstart_url = protocol + "//www.broadinstitute.org/igv/projects/current/igv.php";

    if (queryString) {
        webstart_url += "?" + queryString;
    }

    // determine if webstart is available - code taken from sun site
    var userAgent = navigator.userAgent.toLowerCase();
    // user is running windows
    if (userAgent.indexOf("msie") != -1 && userAgent.indexOf("win") != -1) {
        document.write("<OBJECT " +
            "codeBase=http://java.sun.com/update/1.5.0/jinstall-1_5_0_05-windows-i586.cab " +
            "classid=clsid:5852F5ED-8BF4-11D4-A245-0080C6F74284 height=0 width=0>");
        document.write("<PARAM name=app VALUE=" + webstart_url + ">");
        document.write("<PARAM NAME=back VALUE=true>");
        // alternate html for browsers which cannot instantiate the object
        document.write("<A href=\"http://java.sun.com/j2se/1.5.0/download.html\">Download Java WebStart</A>");
        document.write("</OBJECT>");
    }
    // user is not running windows
    else if (navigator.mimeTypes['application/x-java-jnlp-file']) {
        window.location = webstart_url;
    }
    // user does not have jre installed or lacks appropriate version - direct them to sun download site
    else {
        openWindow("http://jdl.sun.com/webapps/getjava/BrowserRedirect?locale=en&host=java.com",
            "needdownload");
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

getQueryValue = function (name, queryString) {

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