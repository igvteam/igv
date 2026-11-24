# figwig reader: required API

This is the contract between the `figwig` package (the reader of figwig `.fw` files, format version 1)
and an application that draws from it. It lists what IGV requires. An implementation of the package
may offer more, but must provide everything here with the behaviour stated, and IGV must use nothing
else.

The contract is two types: `figwig.FigWigFile`, which the package provides, and `figwig.FigWigByteSource`,
which the application may implement to supply a file's bytes.

## Conventions

These hold for every method unless it says otherwise.

- **Coordinates** are 0-based and half-open: `[start, end)` covers bases `start` to `end - 1`.
- **Track** here always means a track of the figwig file: one signal over the genome, stored in the
  file under a name. It is not an application's track object (IGV's `Track`); how those relate to the
  file's tracks is the application's business.
- **`trackIndex`** identifies a track of the file: `0` to `nTracks() - 1`, in the file's order.
  **`chrom`** identifies a chromosome the same way: `0` to `nChroms() - 1`. Both are fixed for the
  life of the `FigWigFile`.
- **Values** are `double`. A base holds one value per track. A value of zero means the track has
  nothing at that base, and zero is never reported; `NaN` is a value and is reported.
- **Arguments out of range.** A `trackIndex`, `chrom` or `level` outside its range has undefined
  behaviour (in practice an unchecked exception). Coordinate arguments are handled as each method states.
- **Threads.** Every method of `FigWigFile` may be called from any number of threads at once, on the
  same file and for the same track.
- **Errors.** A method fails in one of three ways, and a caller must handle all three:
  - `java.io.IOException`, where the method declares it: the bytes could not be read, the file ended
    early, or the file is corrupt.
  - `java.io.UncheckedIOException`, from methods that do not declare `IOException`: the same causes.
  - `java.lang.IllegalStateException`, from any method that reads the file: the file breaks a rule
    of the format (a checksum, a count, an offset).

  No method throws any other exception because of the file's contents, and none hangs on a corrupt
  or truncated file.

## `FigWigByteSource`

```java
public interface FigWigByteSource extends java.io.Closeable {
    void read(long offset, byte[] out) throws IOException;
    long length() throws IOException;
    boolean isRemote();
}
```

The application implements this for any file the reader should not open itself: a remote file, or one
that needs the application's authentication, signing or redirects.

| Method | Contract for the implementer |
|---|---|
| `read(offset, out)` | Fill all of `out` with the bytes at `offset` to `offset + out.length - 1`. If the file ends first, throw `IOException`; never return with `out` partly filled. Must be safe to call from several threads at once. |
| `length()` | The file's length in bytes, or `0` if it is not known. Called once, while the file is opened. |
| `isRemote()` | `true` if each `read` is a round trip to a server. Must return the same value every time. When `true`, the reader asks for longer stretches per read and issues the reads of one request concurrently. |
| `close()` | Release the source. The reader calls it when the `FigWigFile` is closed, or from `open` if the file cannot be opened. It may be called more than once; calls after the first must do nothing. The application must not close a source it has passed to `open`. |

The reader may ask for up to 2<sup>31</sup> - 1 bytes in one `read`. A `read` that throws fails the call
on `FigWigFile` that needed it; a later call may ask for the same bytes again.

## `FigWigFile`

### Opening and closing

```java
public static FigWigFile open(String path) throws IOException
public static FigWigFile open(String path, FigWigByteSource source) throws IOException
public String path()
public void close() throws IOException
public boolean isClosed()
```

- **`open(path)`** opens a local file. `path` is a file system path, not a URL.
- **`open(path, source)`** opens a file whose bytes come from `source`. `path` is used only as the
  file's name in messages and as the value of `path()`; the reader does not open it.
- Both read and check the file's header and tables before returning, and fail as under *Errors* if the
  file is not a figwig file of version 1. The file owns `source` from the call on: if
  `open(path, source)` fails, it has closed `source` before failing, and no file is left open.
- **`path()`** returns the `path` given to `open`.
- **`close()`** closes the file and its `FigWigByteSource`, whatever is in use (see below). It may be
  called while other threads are calling the file. A call on the file that is running when it is
  closed, or made after, either returns what it would have returned or fails as under *Errors*; it
  never returns wrong values.
- Once the file is closed, a call that needs bytes from the file fails with `IOException` (or
  `UncheckedIOException`) instead of reading them. A `read` on the source that was under way at the
  moment of closing may still complete; its bytes are discarded.
- **`isClosed()`** is `true` once the file has been closed, by `close()` or by the last `release`.

### Tracks in use

```java
public void retain(int trackIndex)
public void release(int trackIndex) throws IOException
```

One figwig file holds many tracks, and an application opens it once and shares the `FigWigFile`
among the objects that draw from it, typically one per track shown. These two methods let each of
those objects say which track of the file it draws, and for how long. A track of the file is
*in use* from a `retain(trackIndex)` until the matching `release(trackIndex)`. The file counts uses: a
track retained twice is in use until it is released twice.

They serve two purposes.

1. **Closing the file.** No one object knows when the file is no longer needed, so the file counts
   its uses and closes itself when the last one ends.
2. **Sharing work.** Tracks of a file are stored together, so what is read from the file for one track
   often serves others. Knowing which tracks are in use lets the reader prepare for the calls that will
   come (`read`, `levelMeans` and `prefetch` with that `trackIndex`) and combine work across them.

- **`retain(trackIndex)`** marks one more use of the track. It returns at once and does not fail. It may
  start reading the file in the background; an error there is not reported until a later call needs
  the same bytes.
- **`release(trackIndex)`** ends one use of the track. If no use of any track of the file remains, it
  closes the file exactly as `close()` does, and fails as `close()` may; otherwise it does nothing
  else and does not fail.
- Every `retain(trackIndex)` must be followed, in time, by exactly one `release(trackIndex)` with
  the same `trackIndex`. A `release` without a `retain` before it has undefined behaviour.
- A file on which `retain` has never been called is closed only by `close()`.
- Using a track without retaining it is allowed and returns the same values. Whether a track is in
  use never changes a result, only how much work the reader does and when.
- A track in use may be used from any number of threads, and by callers that did not retain it.

### File metadata

```java
public int nChroms()
public String chromName(int chrom)
public long chromLength(int chrom)
public int nTracks()
public String trackName(int trackIndex)
public java.util.Map<String, String> trackAttributes(int trackIndex)
public double trackMin(int trackIndex)
public double trackMax(int trackIndex)
```

None of these read the file or fail; the values were read by `open`.

- **`nChroms()`**, **`nTracks()`**: at least 1.
- **`chromName(chrom)`**: the chromosome's name exactly as the file spells it. The reader applies no
  aliases; matching names to a genome is the caller's work.
- **`chromLength(chrom)`**: the chromosome's length in bases.
- **`trackName(trackIndex)`**: the track's `name` from the file's metadata. The reader does not require
  names to be unique.
- **`trackAttributes(trackIndex)`**: every other key of the track's metadata record except `dtype` and
  `group`, each value converted to a string. Never null; empty if there are none. The map is
  read-only, and its iteration order is unspecified.
- **`trackMin(trackIndex)`**, **`trackMax(trackIndex)`**: the least and greatest value of the track as the file's
  statistics table records them. The reader does not compute or check them.

### Values

```java
public Runs read(int trackIndex, int chrom, int start, int end) throws IOException

public static class Runs {
    public int n;
    public int[] starts;
    public int[] ends;
    public double[] values;
}
```

**`read`** returns the track's values over `[start, end)` of the chromosome, exactly as stored.

- `start` below 0 is taken as 0, and `end` past the chromosome's length as the length. If nothing is
  left (`end <= start` after that), the result has `n == 0`.
- The result holds `n` runs. Run `i` says that every base of `[starts[i], ends[i])` has the value
  `values[i]`. Only entries `0` to `n - 1` of the arrays mean anything; the arrays may be longer.
- For every run, `starts[i] < ends[i]`, the run lies within the range asked for (after the clamping
  above), and `values[i]` is not zero. `values[i]` may be `NaN`.
- Runs are in ascending order and do not overlap: `ends[i] <= starts[i + 1]`.
- Every base of the range that holds a nonzero value is in exactly one run; every base in no run
  holds zero.
- Neighbouring bases with the same value are normally reported as one run, but a caller must not
  depend on runs being as long as possible.
- The result belongs to the caller. The reader keeps no reference to it.
- Coordinates are `int`, so `read` reaches only the first 2<sup>31</sup> - 1 bases of a chromosome.

### Summaries

```java
public int summaryBases()
public int levels()
public long levelBases(int level)
public boolean hasSummary(int trackIndex)
public Runs levelMeans(int trackIndex, int level, int chrom, int start, int end)
```

A file may hold its tracks' values summed over bins of fixed size, at several sizes ("levels"). Bins
of a level are laid end to end from base 0 of each chromosome; the last bin of a chromosome is cut
off at the chromosome's length.

- **`levels()`**: the number of levels of the file, at least 1. Levels are numbered from 0. The
  number is the same for every track, including tracks for which `hasSummary` is `false`. May read
  the file the first time it is called.
- **`levelBases(level)`**: the bin size of a level, in bases. `levelBases(0)` is the finest. The
  reader does not check that sizes increase with the level number, so a caller choosing a level must
  compare sizes, not numbers. May read the file the first time it is called.
- **`summaryBases()`**: equal to `levelBases(0)`, but never reads the file or fails.
- **`hasSummary(trackIndex)`**: whether the file holds summaries for this track. A file may hold them for
  some tracks and not others. May read the file the first time it is called for a track.
- **`levelMeans(trackIndex, level, chrom, start, end)`**: the track's mean per bin of the level, for the
  bins that overlap `[start, end)`.
  - Requires `0 <= start < end`. `end` past the chromosome's length is taken as the length.
  - Returns `null` if `hasSummary(trackIndex)` is `false`.
  - Otherwise returns `Runs` with the properties given under *Values*, except these:
    - Runs are not cut to `[start, end)`. Each covers whole bins, so the first may begin before
      `start` and the last end after `end`. No run passes the chromosome's length.
    - `values[i]` is the mean of the track over the bin: the sum of its values at the bin's bases
      divided by the number of bases of the bin within the chromosome. Bases holding zero count
      as zero.
    - A bin whose mean is zero is not reported. Neighbouring bins with the same mean may be
      reported as one run, so a run may span several bins.
  - Fails with `IllegalStateException` if the levels stored for the track are not the file's levels.

### Hint

```java
public void prefetch(int trackIndex, int chrom, int start, int end) throws IOException
```

**`prefetch`** tells the reader that the caller is about to `read` the track over `[start, end)` in
several calls covering parts of it. It changes no result. The reader decides whether to act on it
(by default only when the file's source is remote) and may return at once. It clamps `start` and
`end` as `read` does. A caller need never call it; calling it is correct for any range.

## What is not part of the contract

- **System properties.** The reader reads properties named `figwig.*` to tune itself. Their names,
  defaults and effects may change.
- **Threads.** The reader starts daemon threads of its own. How many, and what they are called, may change.
- **Caching.** What the reader keeps in memory, and for how long, may change. No method here clears it.
- **Other public members.** `FigWigFile` has further public members that IGV's tests use
  (`chromId(String)`, `readDense(...)`) or that report on the reader (`chunkBytesRead()`,
  `chunkReads()`, `clearCache()`, `retain()` and `release()` without a `trackIndex`, `MAGIC`). IGV's own code must not
  depend on them.
- **Everything else in the package** (`FigWigChunks`, `FigWigPool`, `FigWigSummaryIndex`) is internal.
