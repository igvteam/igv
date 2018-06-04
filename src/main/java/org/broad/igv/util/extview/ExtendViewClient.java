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

package org.broad.igv.util.extview;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.Alignment;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;

import java.io.IOException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

/**
 * Port of perl script blatPlot.pl   http://genomewiki.cse.ucsc.edu/index.php/Blat_Scripts
 *
 * @author jrobinso
 *         Date: 11/21/12
 *         Time: 8:28 AM
 */
/** modified from BlatClient.java, jc 20160708 **/

public class ExtendViewClient {

    static int sleepTime = 15 * 1000;  //	#	milli seconds to wait between requests

    public static void postExtendView(final String fName, final int start, final int end, 
                                   final String r_chr, final int r_start, final int r_end) {

        MessageUtils.showMessage("fName :" + fName +" start:" + start + " end:" + end + "<br>" + 
                                 "chr :" + r_chr +" start:" + r_start + " end:" + r_end);
        Map<String, String> params = new HashMap();
        params.put("name", fName); 
        params.put("start", String.valueOf(start));
        params.put("end", String.valueOf(end));
        params.put("chr", r_chr);
        params.put("r_start", String.valueOf(r_start));
        params.put("r_end", String.valueOf(r_end));
        
        String $url = PreferencesManager.getPreferences().get(Constants.EXTVIEW_URL);
        String urlString = ($url + "/FeatureRange/");
        try {
            String result = HttpUtils.getInstance().doPost(HttpUtils.createURL(urlString), params);
            MessageUtils.showMessage("results:" + result);
        } catch (IOException e1) {
            MessageUtils.showErrorMessage("Error in opening extend view: FeatureRange ", e1);
        }

    }

    public static void postExtendView(Alignment aln) {
        
        String chr = aln.getChr();
        int start = aln.getAlignmentStart();
        int end = aln.getAlignmentEnd();

        start -= 10000;
        end += 10000;
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        final byte[] seqBytes = genome.getSequence(chr, start, end);
        String userSeq = new String(seqBytes);

        /*        
        MessageUtils.showMessage("<html>doSVisQuery called with aln" + 
                                 "<br> chr: "+ aln.getChr() + 
                                 "<br> start: "+ aln.getAlignmentStart() + 
                                 "<br> end: " + aln.getAlignmentEnd() +
                                 "<br> strand: " + aln.getReadStrand() +
                                 "<br> readName: " + aln.getReadName() +
                                 "</html>");
        */

        Map<String, String> params = new HashMap();

        params.put("chr", chr);
        params.put("start", String.valueOf(start));
        params.put("end", String.valueOf(end));
        params.put("strand", String.valueOf(aln.getReadStrand()));
        params.put("read_name", aln.getReadName());
        params.put("read_seq", aln.getReadSequence());
        params.put("ref_seq", userSeq);
    

        String $url = PreferencesManager.getPreferences().get(Constants.EXTVIEW_URL);
        String urlString = ($url + "/ExamineReadAlignment/");
        try {
            String result = HttpUtils.getInstance().doPost(HttpUtils.createURL(urlString), params);
            //MessageUtils.showMessage("results:" + result);
        } catch (IOException e1) {
            MessageUtils.showErrorMessage("Error in opening extend view: ExamineReadAlignment", e1);
        }
    }
}

