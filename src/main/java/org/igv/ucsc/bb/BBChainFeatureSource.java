package org.igv.ucsc.bb;

import org.igv.feature.Exon;
import org.igv.feature.IGVFeature;
import org.igv.feature.PSLRecord;
import org.igv.feature.genome.Genome;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.CompletionException;

/**
 * Feature source for UCSC "bigChain" tracks.  Chains are read from the bigChain file, and their aligned blocks from
 * the companion "bigLink" file, where each record is one ungapped block named by its chain id.  Only the links
 * overlapping the queried interval are attached, as in the UCSC browser, so blocks of a chain extending beyond the
 * interval are not loaded.  Where the interval cuts a chain inside a gap, a 1 bp placeholder block is added at that
 * interval edge, as UCSC does, so the gap is drawn as a line (the renderer only connects consecutive blocks, and
 * draws a chain with no blocks at all as a solid block).
 * The chain and link queries are independent and are issued concurrently, as both are typically remote.
 */
public class BBChainFeatureSource extends BBFeatureSource {

    private final BBFeatureSource linkSource;    // Optional, can be null

    public BBChainFeatureSource(BBFile chainReader, BBFile linkReader, Genome genome) throws IOException {
        super(chainReader, genome);
        this.linkSource = linkReader == null ? null : new BBFeatureSource(linkReader, genome);
    }

    @Override
    public Iterator<IGVFeature> getFeatures(String chr, int start, int end) throws IOException {

        CompletableFuture<Iterator<IGVFeature>> linkFuture = linkSource == null ? null :
                CompletableFuture.supplyAsync(() -> {
                    try {
                        return linkSource.getFeatures(chr, start, end);
                    } catch (IOException e) {
                        throw new CompletionException(e);
                    }
                });

        List<IGVFeature> chains = new ArrayList<>();
        Map<String, PSLRecord> chainMap = new HashMap<>();
        Iterator<IGVFeature> iter = super.getFeatures(chr, start, end);
        while (iter.hasNext()) {
            IGVFeature f = iter.next();
            chains.add(f);
            if (f instanceof PSLRecord) {
                chainMap.put(f.getIdentifier(), (PSLRecord) f);
            }
        }

        if (linkFuture != null) {
            Iterator<IGVFeature> linkIter;
            try {
                linkIter = linkFuture.join();
            } catch (CompletionException e) {
                if (e.getCause() instanceof IOException) {
                    throw (IOException) e.getCause();
                }
                throw e;
            }
            while (linkIter.hasNext()) {
                IGVFeature link = linkIter.next();
                PSLRecord chain = chainMap.get(link.getName());
                if (chain != null) {
                    chain.addExon(new Exon(chr, link.getStart(), link.getEnd(), chain.getStrand()));
                }
            }
            for (PSLRecord chain : chainMap.values()) {
                int s = Math.max(chain.getStart(), start);
                int e = Math.min(chain.getEnd(), end);
                List<Exon> exons = chain.getExons();
                if (exons == null || exons.isEmpty()) {
                    chain.addExon(new Exon(chr, s, s + 1, chain.getStrand()));
                    exons = chain.getExons();
                } else if (exons.get(0).getStart() > s) {
                    exons.add(0, new Exon(chr, s, s + 1, chain.getStrand()));   // Links are in position order
                }
                if (exons.get(exons.size() - 1).getEnd() < e) {
                    chain.addExon(new Exon(chr, e - 1, e, chain.getStrand()));
                }
            }
        }

        return chains.iterator();
    }

    /**
     * The link file is typically much denser than the chain file, so size the visibility window by both.
     */
    @Override
    public int getFeatureWindowSize() {
        int chainWindow = super.getFeatureWindowSize();
        int linkWindow = linkSource == null ? -1 : linkSource.getFeatureWindowSize();
        return chainWindow < 0 ? linkWindow : linkWindow < 0 ? chainWindow : Math.min(chainWindow, linkWindow);
    }

    @Override
    public void close() {
        super.close();
        if (linkSource != null) {
            linkSource.close();
        }
    }
}
