package org.broad.igv.util.ucsc;

import org.broad.igv.util.HttpUtils;

import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;

/**
 *
 * Port of perl script blatPlot.pl   http://genomewiki.cse.ucsc.edu/index.php/Blat_Scripts
 * @author jrobinso
 *         Date: 11/21/12
 *         Time: 8:28 AM
 */
public class BlatClient {

    static int  $sleepTime = 120;  //	#	seconds to sleep after request returned

    static void Usage() {
        System.out.println("usage: BlatBot <organism> <db> <searchType> <sortOrder>");
        System.out.println( " <outputType> <querySequence>");
        System.out.println( "\tSpecify organism using the common name with first letter");
        System.out.println( "capitalized.");
        System.out.println( "\te.g. Human, Mouse, Rat etc.");
        System.out.println( "\tDb is database or assembly name e.g hg17, mm5, rn3 etc.");
        System.out.println( "\tsearchType can be BLATGuess, DNA, RNA, transDNA or transRNA");
        System.out.println( "\tsortOrder can be query,score; query,start; chrom,score");
        System.out.println( "\tchrom,start; score.");
        System.out.println( "\toutputType can be pslNoHeader, psl or hyperlink.");
        System.out.println( "\tblats will be run in groups of $batchCount sequences, all");
    }


    public static void main(String [] args) throws IOException {

    if (args.length != 6) {
        Usage();
        System.exit(255);
    }

        String org = args[0];
        String db = args[1];
        String searchType = args[2];
        String sortOrder = args[3];
        String outputType = args[4];
        String userSeq = args[5];


        if (searchType.equals("BLATGuess")) {
        searchType = "Blat's Guess";
    } else if (searchType .equals("transDNA")) {
        searchType = "translated DNA";
    } else if (searchType .equals("transRNA")) {
        searchType = "translated RNA";
    } else if (searchType .equals( "DNA") || (searchType.equals("RNA"))) {
    } else {
        System.out.println("ERROR: have not specified an acceptable search type - it should be BLATGuess, transDNA, transRNA, DNA or RNA.");
        Usage();
        System.exit(255);
    }
    if (outputType.equals("pslNoHeader")) {
        outputType = "psl no header";
    } else if (outputType.equals( "psl") || outputType.equals("hyperlink")) {
    } else {
        System.out.println("ERROR: have not specified an acceptable output type - it should be pslNoHeader, psl or hyperlink.");
        Usage();
        System.exit(255);
    }




    //$response;
    String  $url = "http://genome.cse.ucsc.edu/cgi-bin/hgBlat";

        //if an hgsid was obtained from the output of the first batch
        //then use this.

        URL url = new URL($url + "?org=" + org + "&db=" + db + "&type=" + searchType + "&sort=" + sortOrder +
                "&output=" + outputType + "&userSeq=" + userSeq); // + "&hgsid=" + hgsid);

        String result =   HttpUtils.getInstance().getContentsAsString(url);

            System.out.println(result);

    }
}

//
//sub getHgsid {
//        my ($file, $line, $hgsid, $found);
//$file = shift;
//$hgsid = 0;
//$found = 0;
//# open output file and parse out hgsid and reuse this for each batch
//        # to stop the sessionDb filling up
//        print STDERR "Opening output file $file to retrieve hgsid\n";
//open(FILE, "<$file") || die "Can not open $file: $!\n";
//while (<FILE>) {
//        $line = $_;
//if ($line =~ m/hgsid=([0-9]+)/ && (! $found)) {
//        $hgsid = $1;
//$found = 1;
//}
//        if ($found) {
//        close FILE;
//return $hgsid;
//}
//        }
//        return $hgsid;
//}
