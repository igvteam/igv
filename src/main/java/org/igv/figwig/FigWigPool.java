package org.igv.figwig;

import java.io.IOException;
import java.io.InterruptedIOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;

/**
 * The threads every open figwig file shares: a pool that decodes the parts of a read alongside the
 * thread that asked for it, and over HTTP a pool for requests sent at once.
 */
final class FigWigPool {

    private FigWigPool() {
    }

    // Threads that load and decode the chunks of one read, the reading thread included: a view of one
    // track then uses more than one core.
    static final int DECODE_THREADS = Integer.getInteger("figwig.decodeThreads",
            Math.max(1, Math.min(16, Runtime.getRuntime().availableProcessors())));
    static final ExecutorService DECODER = DECODE_THREADS <= 1 ? null :
            Executors.newFixedThreadPool(DECODE_THREADS - 1, r -> {
                Thread t = new Thread(r, "figwig-decode");
                t.setDaemon(true);
                return t;
            });
    /** A task the pool or the thread that made it runs, whichever takes it first. */
    static final class Shared<T> implements Runnable {
        final Callable<T> task;
        final AtomicBoolean taken = new AtomicBoolean();
        final CountDownLatch done = new CountDownLatch(1);
        T result;
        Throwable error;

        Shared(Callable<T> task) {
            this.task = task;
        }

        public void run() {
            if (!taken.compareAndSet(false, true)) {
                return;
            }
            try {
                result = task.call();
            } catch (Throwable t) {
                error = t;
            } finally {
                done.countDown();
            }
        }
    }
    /**
     * Run the tasks on the decoding pool and this thread; their results in order.  This thread runs
     * every task the pool has not started, so it never waits on a queued task, and a task started on
     * the pool never waits on one either: no part holds a chunk it has claimed while waiting.
     */
    static <T> List<T> parallel(List<Callable<T>> tasks) throws IOException {
        return parallel(tasks, DECODER);
    }
    static <T> List<T> parallel(List<Callable<T>> tasks,
                                ExecutorService pool) throws IOException {
        List<Shared<T>> shared = new ArrayList<>(tasks.size());
        for (Callable<T> t : tasks) {
            shared.add(new Shared<>(t));
        }
        for (int i = 1; i < shared.size(); i++) {
            pool.execute(shared.get(i));
        }
        for (Shared<T> s : shared) {
            s.run();
        }
        List<T> out = new ArrayList<>(shared.size());
        for (Shared<T> s : shared) {
            try {
                s.done.await();
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                throw new InterruptedIOException("figwig: interrupted");
            }
            if (s.error != null) {
                Throwable t = s.error;
                if (t instanceof IOException) throw (IOException) t;
                if (t instanceof RuntimeException) throw (RuntimeException) t;
                if (t instanceof Error) throw (Error) t;
                throw new IOException(t);
            }
            out.add(s.result);
        }
        return out;
    }
    // Over HTTP, parts of loads and background loads wait on the network, not the CPU: their own pool,
    // large enough that every request of a view goes out at once.
    static final ExecutorService PREFETCHER = Executors.newFixedThreadPool(
            Integer.getInteger("figwig.httpThreads", 32), r -> {
        Thread t = new Thread(r, "figwig-http");
        t.setDaemon(true);
        return t;
    });
}
