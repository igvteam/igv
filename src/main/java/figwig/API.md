# figwig reader: required API

This is the contract between the `figwig` package (the reader of figwig `.fw` files, format version 1)
and an application that draws from it. It lists what IGV requires. An implementation of the package
may offer more, but must provide everything here with the behaviour stated, and IGV must use nothing
else.

The contract is two types: `figwig.FigWigFile`, which the package provides, and `figwig.FigWigSource`,
which the application may implement to supply a file's bytes.

## Conventions

These hold for every method unless it says otherwise.

- **Coordinates** are 0-based and half-open: `[start, end)` covers bases `start` to `end - 1`.
- **Tracks** are numbered `0` to `nTracks() - 1`, and **chromosomes** `0` to `nChroms() - 1`, both in
  the file's order. The numbers are fixed for the life of the `FigWigFile`.
- **Values** are `double`. A base holds one value per track. A value of zero means the track has
  nothing at that base, and zero is never reported; `NaN` is a value and is reported.
- **Arguments out of range.** A track, chromosome or level index outside its range has undefined
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

## `FigWigSource`

```java
public interface FigWigSource extends java.io.Closeable {
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
| `close()` | Release the source. The reader calls it once, from `FigWigFile.close()`. |

The reader may ask for up to 2<sup>31</sup> - 1 bytes in one `read`. A `read` that throws fails the call
on `FigWigFile` that needed it; a later call may ask for the same bytes again.

## `FigWigFile`

### Opening and closing

```java
public static FigWigFile open(String path) throws IOException
public static FigWigFile open(String path, FigWigSource source) throws IOException
public String path()
public void retain(int track)
public void release(int track) throws IOException
public void close() throws IOException
public boolean isClosed()
```

- **`open(path)`** opens a local file. `path` is a file system path, not a URL.
- **`open(path, source)`** opens a file whose bytes come from `source`. `path` is used only as the
  file's name in messages and as the value of `path()`; the reader does not open it.
- Both read and check the file's header and tables before returning, and fail as under *Errors* if the
  file is not a figwig file of version 1. If `open(path, source)` fails, `source` has not been closed
  and still belongs to the caller.
- **`path()`** returns the `path` given to `open`.
- **`retain(track)`** declares that the caller will read `track`, and **`release(track)`** that it has
  finished. Each `retain(track)` must be matched by exactly one `release(track)`. The `release` that
  brings the number of retained tracks of the file to zero closes the file, as `close()` does.
  `retain` returns at once and never fails; it may start reading the file in the background.
- A caller that reads without retaining gets the same values. Retaining is what lets the reader
  share work between the tracks in use, and is the only way a file is closed without `close()`.
- **`close()`** closes the file whatever has been retained, and closes its `FigWigSource`. It may be
  called while other threads are reading. A call on the file that is running when it is closed, or
  made after, either returns what it would have returned or fails as under *Errors*; it never
  returns wrong values.
- **`isClosed()`** is `true` once `close()` has been called, directly or by the last `release`.

### File metadata

```java
public int nChroms()
public String chromName(int chrom)
public long chromLength(int chrom)
public int nTracks()
public String trackName(int track)
public java.util.Map<String, String> trackAttributes(int track)
public double trackMin(int track)
public double trackMax(int track)
```

None of these read the file or fail; the values were read by `open`.

- **`nChroms()`**, **`nTracks()`**: at least 1.
- **`chromName(chrom)`**: the chromosome's name exactly as the file spells it. The reader applies no
  aliases; matching names to a genome is the caller's work.
- **`chromLength(chrom)`**: the chromosome's length in bases.
- **`trackName(track)`**: the track's `name` from the file's metadata. The reader does not require
  names to be unique.
- **`trackAttributes(track)`**: every other key of the track's metadata record except `dtype` and
  `group`, each value converted to a string. Never null; empty if there are none. The map is
  read-only, and its iteration order is unspecified.
- **`trackMin(track)`**, **`trackMax(track)`**: the least and greatest value of the track as the file's
  statistics table records them. The reader does not compute or check them.

### Values

```java
public Runs read(int track, int chrom, int start, int end) throws IOException

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
public boolean hasSummary(int track)
public Runs levelMeans(int track, int level, int chrom, int start, int end)
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
- **`hasSummary(track)`**: whether the file holds summaries for this track. A file may hold them for
  some tracks and not others. May read the file the first time it is called for a track.
- **`levelMeans(track, level, chrom, start, end)`**: the track's mean per bin of the level, for the
  bins that overlap `[start, end)`.
  - Requires `0 <= start < end`. `end` past the chromosome's length is taken as the length.
  - Returns `null` if `hasSummary(track)` is `false`.
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
public void prefetch(int track, int chrom, int start, int end) throws IOException
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
  `chunkReads()`, `clearCache()`, `retain()`, `release()`, `MAGIC`). IGV's own code must not
  depend on them.
- **Everything else in the package** (`FigWigChunks`, `FigWigPool`, `FigWigSummaryIndex`) is internal.
