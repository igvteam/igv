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
