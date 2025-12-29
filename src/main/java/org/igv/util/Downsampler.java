package org.igv.util;

import java.util.Arrays;
import java.util.Random;

/**
 * Downsample input list of objects using reservoir sampling.
 */
public class Downsampler<T> {

    private final Random RAND = new Random();

    public T[] sample(T[] input, int max) {
        if (input.length < max) {
            return input;
        } else {
            T[] downsampled = Arrays.copyOf(input, max);
            for (int i = max; i < input.length; i++) {
                double samplingProb = ((double) max) / (i + 1);
                if (RAND.nextDouble() < samplingProb) {
                    int idx = (int) (RAND.nextDouble() * (max - 1));
                    downsampled[idx] = input[i];
                }
            }
            return downsampled;
        }
    }

}
