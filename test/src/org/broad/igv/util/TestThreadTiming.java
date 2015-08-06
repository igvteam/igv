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

/*
 * Simple class to time the overhead of splitting a long running job into
 * multiple threads.  On a Macbook Pro the overhead is ~ 20%.  
 */

package org.broad.igv.util;

import org.junit.Ignore;

/**
 * @author jrobinso
 */
@Ignore
public class TestThreadTiming {

    static class Compute implements Runnable {

        int nRepeats;

        Compute(int nRepeats) {
            this.nRepeats = nRepeats;
        }

        public void run() {
            double sum = 0.0;
            for (int i = 0; i < nRepeats; i++) {
                //sum += i * Math.random();
                sum += 1000000.5d;
            }
            System.out.println("Done: " + sum);
        }

    }

    public static void main(String[] args) {

        int nRepeats = 10000000;

        // First test -- one thread
        Compute c1 = new Compute(nRepeats);
        long t0 = System.currentTimeMillis();
        c1.run();
        long dt = System.currentTimeMillis() - t0;
        System.out.println("Single threaded time: " + dt);

        // Second pass -- 10 threads
        t0 = System.currentTimeMillis();
        int nThreads = 10;
        Thread[] workerThreads = new Thread[10];
        for (int i = 0; i < nThreads; i++) {
            Compute c2 = new Compute(nRepeats / nThreads);
            workerThreads[i] = new Thread(c2);
            workerThreads[i].start();
        }

        // Wait for workers to complete
        try {
            for (int i = 0; i < nThreads; i++) {
                workerThreads[i].join();
            }

        } catch (InterruptedException ex) {
            ex.printStackTrace();
        }

        dt = System.currentTimeMillis() - t0;
        System.out.println("Multi threaded time:" + dt);
    }

}
