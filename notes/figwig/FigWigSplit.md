# Splitting the figwig reader into its own repository

Notes on PR [#1900](https://github.com/igvteam/igv/pull/1900), "Read figwig (.fw) multi-track signal
files" (branch `jmschrei-figwig-phase0`). The question: can the figwig reader and its tests move to a
separate repository that builds a jar IGV depends on?

Status, 2026-11-24: done. The reader is its own project, https://github.com/igvteam/figwig (package
`org.figwig`, MIT), and release v0.1.0 has `figwig-0.1.0.jar` attached. The IGV branch
`figwig-external-jar` (local) depends on that jar and no longer holds the reader. Nothing has been
proposed on the PR yet.

## The library and how IGV gets it

- **Repository:** `igvteam/figwig`, intended to move to the lab that owns the format. It holds the
  reader, `FigWigByteSource`, `FigWigReader` (a command-line reader after the Python API), their
  tests and test data, `docs/API.md` (the contract for what IGV uses) and a README.
- **Build:** a minimal Gradle build whose only dependency is `org.json`. `./gradlew build` makes the
  jar; `./gradlew installDist` adds a `figwig` launcher.
- **Releasing:** publishing a GitHub release tagged `v1.2.3` runs a workflow that tests, builds
  `figwig-1.2.3.jar` and attaches it. The version comes from the tag; no file is edited.
- **IGV's build:** an `ivy` repository in `build.gradle` that maps `igvteam:figwig:0.1.0` to the
  release's download URL. No credentials or third-party service. It carries no dependency metadata,
  which does not matter while the library needs only `org.json`.
- **Not chosen:** GitHub Packages (needs a token even to download public packages), JitPack (a
  third-party service in IGV's build), Maven Central (not ours to do for the lab; the swap later is
  one repository block and one dependency line).
- **Result:** on `figwig-external-jar` the full suite passes with 689 tests (the 15 reader tests now
  run in the library's own repository), and `createDist` puts `figwig-0.1.0.jar` in `lib/`.
- **Test data:** the `.fw` files stay in IGV's `test/data/figwig` for its 13 tests; the expected
  tables moved with the library.

The sections below describe the in-tree simulation on `figwig-split-simulation` that preceded this,
where the reader was a top-level `figwig` package inside IGV.

## Summary

Feasible, and confirmed by doing it. The reader now sits in a top-level `figwig` package that imports
nothing from `org.igv` or htsjdk, and `javac` compiles it alone against the `org.json` jar. The IGV
side is three classes (about 620 lines) plus the hooks in `TrackLoader` and `FileFormatUtils`.

- All 24 figwig tests pass, as they did before the change.
- The full IGV suite passes: 700 tests, 0 failures, 37 skipped.
- The change was 13 files, +385 / -208, most of it accessors and the I/O seam.

Still a package boundary only: there is no separate Gradle module or jar, the test data is not split,
and performance has not been re-measured.

## Where the line falls

| Piece | Lines | Side | Why |
|---|---|---|---|
| `FigWigFile`, `FigWigChunks`, `FigWigPool`, `FigWigSummaryIndex` | 2,624 | Library | The format: depend only on the JDK and `org.json` |
| `FigWigByteSource` (new) | 27 | Library | The interface the application implements to supply bytes |
| `FigWigDataSource` | 464 | IGV | Extends `AbstractDataSource`; uses `Genome`, `LocusScore`, `Globals`, `FrameManager`, `LRUCache` |
| `FigWigTrackChooser` | 81 | IGV | Swing UI on `IGV`, `TrackChooser`, `FileRecord` |
| `IGVFigWigByteSource` (new) | 78 | IGV | Adapts IGV's streams to `FigWigByteSource` |
| `FigWigFileTest` | 465, 11 tests | Library | Reader against the figwig package's expected values; no IGV imports |
| `FigWigDataSourceTest`, `FigWigEdgeCaseTest`, `FigWigRemoteTest` | 555, 13 tests | IGV | Go through `TrackLoader`, a genome, or IGV's HTTP stack |

Two tests moved from `FigWigEdgeCaseTest` to the library's `FigWigFileTest` because they only exercise
the reader: ranges reaching past a chromosome, and reads racing the file's closing. The IGV-side tests
no longer reference the library's test class.

## The I/O seam

The reader reached into IGV in two places: `IGVSeekableStreamFactory.getInstance().getStreamFor(path)`
(two call sites) and `FileUtils.isRemote`. Both are replaced by one interface in the library:

```java
public interface FigWigByteSource extends Closeable {
    void read(long offset, byte[] out) throws IOException;  // called from several threads at once
    long length() throws IOException;                       // 0 when not known
    boolean isRemote();                                     // each read is a round trip
}
```

- `FigWigFile.open(path)` opens a local file itself, with positional `FileChannel` reads as before.
- `FigWigFile.open(path, source)` reads everything through the application's source.
- `IGVFigWigByteSource.open(path)` is what `TrackLoader` calls: a local file goes to the reader directly,
  anything else through IGV's streams. That keeps OAuth, S3/GS signing and redirects on the IGV side.
- For remote paths the adapter opens one stream per read, so loader threads still read concurrently;
  other non-local paths share one stream under a lock. Both behaviours are unchanged from the PR.
- The `figwig.noPositional` switch moved to the IGV adapter, since it chooses between the two routes.

A single thread-safe source turned out simpler than the stream factory first considered: the library
never sees a stream, so it needs no htsjdk types at all.

## The API the IGV side uses

`FigWigFile` has no public fields any more, so the compiler stops IGV code from using anything but
methods. What `FigWigDataSource`, `FigWigTrackChooser` and `TrackLoader` call:

| Kind | Methods |
|---|---|
| Opening | `open(path)`, `open(path, source)` |
| File metadata | `path()`, `nChroms()`, `chromName(c)`, `chromLength(c)`, `nTracks()`, `trackName(t)`, `trackAttributes(t)`, `trackMin(t)`, `trackMax(t)` |
| Per-base data | `read(track, chrom, start, end)` returning `Runs` |
| Summaries | `levels()`, `levelBases(k)`, `levelMeans(...)`, `hasSummary(track)`, `summaryBases()` |
| Hint | `prefetch(track, chrom, start, end)` |
| Lifecycle | `retain(track)`, `release(track)`, `close()` |

The IGV-side tests also use `chromId(name)`, `readDense(...)` and `isClosed()`. Other methods that were
already public were left public (`chunkBytesRead()`, `chunkReads()`, `clearCache()`, the no-argument
`retain()` / `release()`).

Points that needed a decision:

- **`summaryBases()` replaces the `tile` field.** The data source needs that number twice: to skip
  summaries finer than the view (`FigWigDataSource.java:107`), and as its threshold for summarizing
  tiles in parallel (`:267`). It is still the tile size underneath; the method only stops IGV from
  depending on the field. A format with a variable or per-group tile size would need a different answer.
- **`prefetch` stays public, as a hint.** It cannot move inside `read`: the data source prefetches the
  whole range of uncached tiles before reading them one at a time, and `read` cannot know that range.
  What moved into the reader is the policy (remote files only, and the `figwig.prefetch` property), so
  `isHttp()` is no longer public.
- **`Runs` keeps its public arrays** (`n`, `starts`, `ends`, `values`). It is the result struct of the
  data source's hot loops; only its mutators are hidden. This is the one place a field layout is still
  part of the contract.
- **The data source has its own pool.** It used the package-private `FigWigPool.DECODE_THREADS` and
  `FigWigPool.parallel`. It now runs `FutureTask`s on a small `figwig-summary` pool, each run by the
  pool or the requesting thread, whichever takes it first, as `FigWigPool.parallel` does. The pools
  are separate, so summarizing and decoding together can use up to twice the threads as before.

## Why `FigWigDataSource` stays in IGV

It is the IGV half of the adapter, not part of the format reader. Moving it would make the library
depend on `igv.jar` while IGV depends on the library.

- **Built from IGV types.** Extends `AbstractDataSource`, implements `DataSource`; inputs and outputs
  are `Genome`, `LocusScore`, `DataTile`, `BasicScore`, `CompositeScore`, `WindowFunction`.
- **Its logic is IGV's display model.** It maps IGV zoom levels to figwig summary levels using the
  700-bin tile convention, builds the whole-genome view in genome coordinates, resolves chromosome
  aliases through `Genome`, and checks `FrameManager.isGeneListMode()`.
- **About half of it re-implements `AbstractDataSource` internals.** The `summaryTiles` / `summarize` /
  `composite` block reproduces `getSummaryScoresForRange` and `computeSummaryTile` "step for step" for
  speed. If `AbstractDataSource` changes, this copy has to change with it, which is only practical in
  the same repository.

Review point independent of the split: that 230-line parallel copy of core summarisation logic is the
part of the PR that puts the most ongoing maintenance on IGV. The PR tests that it matches
`AbstractDataSource` score for score, which catches drift, but it is worth asking whether the speed
gain justifies the duplication.

## Is the data source isolated from format changes?

Yes, for any revision that keeps the per-base and mean-summary model. Nothing in the API above exposes
bytes, chunks, codecs or offsets, and `FigWigChunks` / `FigWigSummaryIndex` are never referenced from
the IGV side. Changes to chunk kinds, element types, codecs, track ordering, the summary index or the
header layout stay inside the library.

What still couples the two:

- **Semantic assumptions in the caller.** Summaries are means only (so only `mean` and `none` are
  offered), zero means "no data" and is omitted from runs, and values are narrowed from `double` to
  `float`. New summary kinds or a distinct missing-value marker would need data source changes to be
  used, though the existing code would keep working.
- **`summaryBases()` and the `Runs` arrays**, as described above.
- **Tuning properties.** The library reads about 20 `figwig.*` system properties. They are not API in
  the Java sense, but anyone setting them depends on the names.

## Build and packaging

- No `module-info.java`; the launch scripts build the classpath from everything in `lib/`; `createDist`
  copies all Maven-resolved jars. A new dependency is one line in `build.gradle`.
- The real cost is publishing. The jar has to come from somewhere `build.gradle` resolves: Maven
  Central, or the Sonatype snapshots repository already configured. That means a release process,
  versioning, and a two-step workflow for any reader fix (release the library, then bump IGV).
- Test data (`test/data/figwig`, 13 `.fw` files and their expected tables, 691 KB) would be split or
  duplicated. The library's tests need all of it; the IGV-side tests read the `.fw` files but not the
  expected tables, and `FigWigRemoteTest` serves every `.fw` file in the directory.

## Open

1. **Who owns the new repository?** The split is worth doing if someone other than the IGV team owns
   the reader. The format is new and still evolving, and 2,600 lines of dense concurrent decoding with
   its own tuning properties is a lot to adopt into a mature codebase. If the IGV team would maintain
   the library anyway, the release overhead buys little and in-tree is simpler.
2. **Performance.** The PR's timings were taken with the data source sharing the reader's decode pool.
   The separate summary pool should be re-measured before the split is taken as free.
3. **Accessor names.** `chromName(c)`, `trackName(t)` and the rest follow the reader's existing style
   (`nTracks()`, `levels()`). They are the author's to name if the library goes ahead.
4. **Whether to push the simulation branch**, so the PR author can see the diff.

Settled by the simulation: the library does not need htsjdk. The source interface made that the
simpler option as well as the lighter one.
