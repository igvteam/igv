package org.igv.alignment;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import org.junit.Assert;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;

/**
 * Tests of the insert size ("TLEN") threshold estimate.
 */
public class PEStatsTest {

    /**
     * A region dense with discordant pairs -- a large deletion, say -- must not drag the maximum threshold
     * past every read worth flagging.  The trimmed sample, and so the thresholds, should be identical to
     * those of the same library with no discordant population at all.
     */
    @Test
    public void testDiscordantPairsDoNotInflateThresholds() {

        PEStats clean = new PEStats("clean");
        clean.merge(concordantLoad());
        clean.computeInsertSize(0.5, 99.5);

        PEStats contaminated = new PEStats("contaminated");
        PEStats load = concordantLoad();
        for (int i = 0; i < 20; i++) {
            load.update(properPair(50000));
        }
        contaminated.merge(load);
        contaminated.computeInsertSize(0.5, 99.5);

        Assert.assertTrue("maximum reflects the library", clean.getMaxThreshold() < 1000);
        Assert.assertEquals("maximum unaffected by outliers", clean.getMaxThreshold(), contaminated.getMaxThreshold());
        Assert.assertEquals("minimum unaffected by outliers", clean.getMinThreshold(), contaminated.getMinThreshold());
    }

    /**
     * No single load may set the thresholds for the session.  A deep load supplies far more proper pairs
     * than the pool holds, so without a per load cap the region visited first would own the pool and a
     * later region could not move the estimate.
     */
    @Test
    public void testNoSingleLoadDominatesThePool() {

        PEStats stats = new PEStats("lib");
        stats.merge(uniformLoad(50000, 300));    // deep region
        stats.merge(uniformLoad(1000, 5000));    // shallow region, different insert size

        // Each load contributes at most 1000, so the second is half the pool, not 2% of it
        stats.computeInsertSize(0.5, 95);
        Assert.assertEquals(5000, stats.getMaxThreshold());
    }

    /**
     * A minimum percentile of zero means no read is to be colored as a small insert.  It is not a
     * percentile the estimate can be taken at, so it must be honored rather than passed to the estimator.
     */
    @Test
    public void testZeroMinPercentileMeansNoMinimum() {

        PEStats stats = new PEStats("lib");
        stats.merge(concordantLoad());
        stats.computeInsertSize(0, 99.5);

        Assert.assertEquals(0, stats.getMinThreshold());
        Assert.assertTrue("maximum still estimated", stats.getMaxThreshold() > 0);

        // Also before there is a sample to estimate from
        PEStats unsampled = new PEStats("lib");
        unsampled.computeInsertSize(0, 99.5);
        Assert.assertEquals(0, unsampled.getMinThreshold());
    }

    /**
     * 980 proper pairs spread over 250-350 bp, as accumulated by a single load.
     */
    private static PEStats concordantLoad() {
        PEStats load = PEStats.forLoad("lib");
        for (int i = 0; i < 980; i++) {
            load.update(properPair(250 + (i % 101)));
        }
        return load;
    }

    private static PEStats uniformLoad(int count, int insertSize) {
        PEStats load = PEStats.forLoad("lib");
        for (int i = 0; i < count; i++) {
            load.update(properPair(insertSize));
        }
        return load;
    }

    private static final Map<Integer, Alignment> alignmentCache = new HashMap<>();

    /**
     * A properly paired FR alignment with the given insert size.  Cached -- the statistics read only the
     * insert size and pair orientation, so a single instance can stand in for many reads.
     */
    private static Alignment properPair(int insertSize) {

        return alignmentCache.computeIfAbsent(insertSize, isize -> {

            SAMFileHeader header = new SAMFileHeader();
            header.addSequence(new SAMSequenceRecord("chr1", 1000000));

            SAMRecord record = new SAMRecord(header);
            record.setReadName("read");
            record.setReferenceName("chr1");
            record.setAlignmentStart(1000);
            record.setCigarString("50M");
            record.setReadString("A".repeat(50));
            record.setReadPairedFlag(true);
            record.setProperPairFlag(true);
            record.setFirstOfPairFlag(true);
            record.setReadNegativeStrandFlag(false);
            record.setMateReferenceName("chr1");
            record.setMateAlignmentStart(1000 + isize - 50);
            record.setMateNegativeStrandFlag(true);
            record.setInferredInsertSize(isize);

            return new SAMAlignment(record);
        });
    }
}
