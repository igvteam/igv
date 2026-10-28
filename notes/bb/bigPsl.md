# bigPsl support: decode bigPsl bigBed records as PSLRecords

## Context

IGV parses ASCII PSL files into `PSLRecord` (extends `BasicFeature`) via `PSLCodec`, and BLAT results use the same class. UCSC's `bigPsl` format is a bigBed (bed12+13, autoSql table name `bigPsl`) carrying PSL records. Today IGV decodes bigPsl through the generic `BBBedCodec`, so tracks get plain `BasicFeature`s with the PSL-specific columns as string attributes, and no PSL score/description. Goal: bigPsl records become `PSLRecord`s built by the same code as the ASCII path.

Decisions made with Jim:
- bigPsl-only columns (`oSequence`, `oCDS`, `seqType`, `oStrand`) are **dropped**; popup matches ASCII PSL.
- `thickStart`/`thickEnd` are **ignored**. ASCII PSL has no CDS columns and `PSLRecord` sets thick = start/end, so bigPsl must do the same: rendering is identical between ASCII PSL and bigPsl for the same alignment, regardless of whether the bigPsl was built with a CDS file.
- Test against the public UCSC GenArk Python molurus bigPsl over HTTP, plus a local hand-built record test.
- Codec selection is on the exact autoSql table name (`bigPsl`, `bigChain`). No known files rename the table, and
  matching on the field set cannot guarantee uniqueness.
- Malformed records fall back to a plain BED feature with a `log.warn`, capped at 10 per file.
- Records are decoded directly into `PSLRecord` (no shared builder with the ASCII codec). Sharing the builder meant
  parsing every record twice; the ~15 lines of score/description code are duplicated instead.
- bigChain is handled the same way: chains become `PSLRecord`s, with blocks from the companion bigLink file.

## Approach

`BBPslCodec` parses the bigPsl columns straight into a `PSLRecord`, populated the same way `PSLCodec.getPslRecord`
populates it from ASCII PSL (same fields, same score formula, same description). Extra-column positions are resolved
once from the autoSql in the constructor. `PSLCodec`'s signatures are untouched; all BLAT callers untouched.

Key bigPsl semantics (confirmed from kent `pslToBigPsl.c` / `bigPsl.c`):
- Target blocks are always on the + reference strand (converter applied `pslRc` when PSL strand[1] was `-`, recording that in `oStrand`). **Do not reverse anything.** IGV draws exactly what bigPsl stores.
- `strand` column = PSL `strand[0]` (query strand). `oChromStart/oChromEnd/oChromSize/oChromStarts` = PSL `qStart/qEnd/qSize/qStarts` verbatim. `chromSize` = `tSize`. `match/misMatch/repMatch/nCount` verbatim.
- `qNumInsert/qBaseInsert/tNumInsert/tBaseInsert` are **not stored**; compute them from blocks like kent `pslComputeInsertCounts`.
- Protein alignments (`seqType == 2`): bigPsl `blockSizes` are in bases (3x residues); query coords are in residues. Exons stay in bases; divide block size by 3 only for the query-gap arithmetic.

## Changes

- `src/main/java/org/igv/ucsc/bb/codecs/BBRecordCodec.java` (new): base class for codecs decoding a bigBed record
  directly into a specialized feature. Holds the autoSql field-index helpers and the malformed-record fallback
  (plain BED feature, warning logged for the first 10 per file).
- `src/main/java/org/igv/ucsc/bb/codecs/BBPslCodec.java` (new): bigPsl -> `PSLRecord`. Column positions of the
  extra fields are resolved once from the autoSql. Blocks come from blockSizes/chromStarts (absolute, + strand);
  qStarts from oChromStarts; insert counts computed from the blocks (kent `pslComputeInsertCounts`, query gaps in
  residues for seqType 2); score and description as in `PSLCodec`.
- `BBCodecFactory`: autoSql table `bigPsl` -> `BBPslCodec`, `bigChain` -> `BBChainCodec`.
- `PSLRecord.setTSize` added and called from `PSLCodec.getPslRecord` (`getText()` previously wrote 0 for tSize).
- bigChain: see below.

## Tests

1. **`src/test/java/org/igv/ucsc/bb/BBPslCodecTest.java`** (local, deterministic). Embed the bigPsl `.as` text as a constant; `BBCodecFactory.getCodec(as, 12)` must return a `BBPslCodec`. Decode a hand-built `BedData` for the second record of `test/data/psl/humanBlat.psl` (`chr6 84433944 84434116`, two blocks `18,9` at `0,163`, qStarts `92,110`, qSize 119, chromSize 170805979, match 27). Assert: result is `PSLRecord`; `getTGapCount()==1`, `getTGapBases()==145`, `getQGapCount()==0`, `getqSize()==119`, `getTSize()==170805979`; two exons at 84433944 (len 18) and 84434107 (len 9); `getScore()` equals (delta 0.01) `new PSLCodec().decode(asciiLine).getScore()`. This covers the shared builder, insert-count logic, `setTSize`, and factory registration.
2. **Remote smoke test in `BBFeatureSourceTest`** (`testBigNarrowpeak` pattern): open `https://hgdownload.soe.ucsc.edu/hubs/GCF/000/186/305/GCF_000186305.1/bbi/GCF_000186305.1_Python_molurus_bivittatus-5.0.2.bigPsl.bb`, query `NW_006532014.1` (largest scaffold, 1.45 Mb) over `0..Integer.MAX_VALUE`; assert count > 0, every feature `instanceof PSLRecord`, `getqSize() > 0` and `getMatch() <= getqSize()`.

## bigChain

bigChain differs from bigPsl structurally: the chain file (`bigChain`, bed6+6: name=chain id, score, strand, tSize,
qName, qSize, qStart, qEnd, chainScore) has **no blocks**. The aligned blocks live in a companion `bigLink` file
(bed4+1: one record per ungapped block, name = chain id, qStart), referenced by the trackDb `linkDataUrl` property or,
by UCSC convention, `<chain>.link.bb`. The UCSC browser queries the link file over the visible window and joins links
to chains by name; IGV does the same.

- `BBChainCodec` decodes a chain record into a `PSLRecord`: name = qName, identifier = chain id, strand, tSize, qSize,
  qStart, qEnd, score = BED score; description shows chain id, chain score and query locus. No exons.
- `BBChainFeatureSource` extends `BBFeatureSource` and holds a second source for the link file.
  `getFeatures(chr, start, end)` reads the chains, then the links over the same interval, and attaches each link as
  an `Exon` of the chain with the matching id. The visibility window is the smaller of the chain and link windows.
  The two queries are independent and run concurrently (the link query on a `CompletableFuture`), saving one round
  of network latency per load on remote hubs.
- `BBLinkCodec` decodes a bigLink record (table `bigLink`) to a minimal feature with only the name, skipping the
  generic BED parsing and the attribute map the chain source never reads. Blocks of a chain outside the queried interval are not loaded (same as UCSC). A chain
  cut by the interval inside a gap gets a 1 bp placeholder block at that interval edge so the gap draws as a line,
  as UCSC does (the renderer only connects consecutive blocks). Without a link file, chains render as a single span. The popup shows query coordinates of `-` strand
  chains converted to the forward strand, as UCSC does.
- `FeatureType.CHAIN`, set by `BBFile.getFeatureType` for autoSql table `bigChain`; `TrackLoader.loadBBFile` builds
  a `FeatureTrack` over a `BBChainFeatureSource` for it. `openLinkFile` uses `locator.getLinkURL()` or the `.link.bb`
  convention, and logs a warning if the link file cannot be opened.
- Hub plumbing: `TrackConfig.linkURL` (also in its JSON in/out), `ResourceLocator.getLinkURL()`,
  `TrackDbHub` maps `bigchain` -> `bigbed` and copies `linkDataUrl`.
- Sessions: `AbstractTrack.marshalJSON` writes `linkURL`; `JSONSessionReader.createLocator` restores `format`,
  `indexURL`/`indexPath` and `linkURL` onto the locator (index restore was previously missing from JSON sessions).
  Without a `linkURL` the loader falls back to the `.link.bb` convention.
- Visibility window: zoomed out past the window the track draws nothing, as for any feature track. The link density
  sets that window, so on a dense chain track it is smaller than the chain file alone would give. If that proves too
  restrictive, the alternative is to skip the link query above some window size and draw chains as plain spans, but
  that is a visible behavior choice rather than a cleanup.
- Test: `BBChainFeatureSourceTest` on the UCSC example files `test/data/bb/bigChain.bb` + `bigChain.link.bb`
  (4 chains, 11 links on chr22_KI270731v1_random).

## Initial display

Alignment tracks are read row by row, so psl, pslx, bigPsl and bigChain tracks load in EXPANDED mode with an
initial height of 300 (`TrackLoader.ALIGNMENT_FEATURE_HEIGHT`) rather than the collapsed default for feature tracks.
`FeatureType.PSL` is set by `BBFile.getFeatureType` for autoSql table `bigPsl` so `loadBBFile` can recognize it; the
tribble path keys on the locator format. The override runs at load time, so session attributes (applied afterwards)
still win, but a hub's `visibility` setting is overridden.

## Verification

```bash
cd /Users/jrobinso/projects/igv
./gradlew test --tests org.igv.ucsc.bb.BBFeatureSourceTest --tests org.igv.ucsc.bb.BBUtilsTest   # baseline before edits
./gradlew test --tests org.igv.ucsc.bb.BBPslCodecTest --tests org.igv.ucsc.bb.BBFeatureSourceTest  # after
./gradlew test   # full regression
```

Manual: `./gradlew createDist`, launch, load the Python bigPsl URL (or any hub bigPsl track) and confirm blocks render and the popup shows `matches = ...` lines. Load `test/data/bb/bigChain.bb` and confirm chains show blocks from the link file.
