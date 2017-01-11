/*
 *  The MIT License (MIT)
 *
 * Copyright (c) 2016 University of California San Diego
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 *
 */

package org.broad.igv.sam;

/**
 * Created by jrobinso on 1/4/17.
 */
//
//import htsjdk.samtools.util.CloseableIterator;
//import htsjdk.samtools.util.RuntimeEOFException;
//import org.broad.igv.sam.reader.BAMReader;
//import org.broad.igv.util.ResourceLocator;
//
//import java.io.IOException;
//import java.util.List;
//import java.util.concurrent.CompletableFuture;
//import java.util.concurrent.ExecutionException;
//import java.util.stream.Collectors;
//

public class TestCF {

    private static final String BAM_URL_STRING = "http://1000genomes.s3.amazonaws.com/phase3/data/HG01879/exome_alignment/HG01879.mapped.ILLUMINA.bwa.ACB.exome.20120522.bam";

//    public List<Alignment> load(String chr, int start, int end) {
//
//        try {
//
//
//            BAMReader reader = new BAMReader(new ResourceLocator(BAM_URL_STRING), true);
//
//            CloseableIterator<PicardAlignment> iter = reader.query(chr, start, end, false);
//
//            List<Alignment> alignments = iter.stream().collect(Collectors.toList());
//            iter.close();
//
//            return alignments;
//        } catch (IOException e) {
//            e.printStackTrace();
//            throw new RuntimeEOFException(e);
//        }
//    }
//
//    public static void main(String[] args) throws IOException, ExecutionException, InterruptedException {
//
//        final String chr = "Y";
//        final int start = 10000000 - 1;
//        final int end = 10004000;
//
//        TestCF cf = new TestCF();
//
////         List<Alignment> alignments = cf.load();
////
////         System.out.println(alignments.size());
//
//        CompletableFuture f = CompletableFuture.supplyAsync(() -> cf.load(chr, start, end)).thenAccept(System.out::println);
//
//
//        f.get();  // To prevent program exit
//    }

}