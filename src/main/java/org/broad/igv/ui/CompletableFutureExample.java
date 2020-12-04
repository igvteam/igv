package org.broad.igv.ui;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class CompletableFutureExample {

    public static void main(String[] args) {
        go();
    }

    static class Task {

        int n;

        public Task(int n) {
            this.n = n;
        }

        void doit() {
            System.out.println("Start " + n);
            try {
                Thread.sleep(1000);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            System.out.println("End " + n);
        }

    }

    private static void go() {

        List<CompletableFuture> futures = new ArrayList();
        for (int i = 0; i < 5; i++) {
            int finalI = i;
            futures.add(CompletableFuture.runAsync(() -> {
                (new Task(finalI)).doit();
            }, threadExecutor));
        }

        final CompletableFuture[] futureArray = futures.toArray(new CompletableFuture[futures.size()]);
        try {
            CompletableFuture.allOf(futureArray).get();
            System.out.println("Done");
            System.exit(0);
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (ExecutionException e) {
            e.printStackTrace();
        } finally {
        }

    }

    private static final ExecutorService threadExecutor = Executors.newFixedThreadPool(5);
}
