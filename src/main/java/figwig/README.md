# figwig for Java

A Java reader of [figwig](https://github.com/jmschrei/figwig) `.fw` files: the per-base values of many
tracks over one genome, in one file. This page covers `FigWigReader`, which reads windows of a file
into arrays and can be run from the command line. It follows the `FigWigReader` of the figwig Python
package; the differences are listed at the end.

`FigWigReader` is a demonstration, not a replacement for the Python package. It reads; it does not
write or convert files.

## Requirements

- Java 21 or later.
- A local `.fw` file of format version 1. Remote files (http, s3, gs) are not read: a URL is taken for
  a file name and fails.
- Two jars on the class path: the one holding the `figwig` package, and `org.json`.

## Building

There is no separate figwig jar yet. The package is built as part of IGV:

```bash
./gradlew createDist
```

This writes `build/IGV-dist/lib/`, which holds `igv.jar` (with the `figwig` package in it) and
`json-20250517.jar`. Those two are all `FigWigReader` needs.

## Running from the command line

```bash
java -cp "build/IGV-dist/lib/igv.jar:build/IGV-dist/lib/json-20250517.jar" figwig.FigWigReader info FILE
java -cp "build/IGV-dist/lib/igv.jar:build/IGV-dist/lib/json-20250517.jar" figwig.FigWigReader read FILE CHROM START WIDTH [TRACK ...]
```

On Windows, separate the two jars with `;` instead of `:`. The examples below shorten the class path
to `$CP`.

### `info`: what a file holds

```
$ java -cp "$CP" figwig.FigWigReader info test/data/figwig/positions_zlib.fw
genome	toy
tile	64
zooms	[4, 16]
chromosomes	3
track	sum	nonzero	min	max	record
0	54.0	54	1.0	1.0	{name=ones, dtype=uint8}
1	63087.0	462	1.0	255.0	{name=random, dtype=uint8}
...
11	-2.125033537799027	156	-3.8994216918945312	2.108288526535034	{name=signal, dtype=float32, assay=ATAC, biosample=K562}
```

- `genome` is the reference the file names, `tile` the bases per tile, and `zooms` the zoom levels the
  file keeps, in tiles per bin.
- Each track line gives the track's index, then its sum, count of nonzero bases, least and greatest
  value as the file records them, then its name, the type its values were written as, and its
  attributes.

### `read`: the values of a window

```
$ java -cp "$CP" figwig.FigWigReader read test/data/figwig/positions_zlib.fw chr1 0 4 signal random
position	signal	random
0	-0.2991911470890045	217.0
1	0.0	163.0
2	0.0	130.0
3	2.108288526535034	69.0
```

- `START` is 0-based, as in BED files: the window covers bases `START` to `START + WIDTH - 1`.
- The output is tab-separated: a header line, then one line per base, with a column per track in the
  order the tracks were named. With no track named, every track of the file is printed.
- Tracks are named by name, exactly as `info` prints them. Quote names the shell would interpret,
  such as `'procap.+'`.
- A base that was never written prints as `0.0`. A base past the end of the chromosome prints as
  `NaN`, and so does a base whose stored value is NaN.
- Every value prints as a decimal number, so a count of 217 prints as `217.0`.

### Things to know

- **Errors** print one line saying what is wrong and the command exits with status 1: an unknown
  chromosome or track name, a start below 0, a file that does not exist, or a file that is not a
  figwig file. For example, `figwig: k562.fw has no chromosome 1`. Wrong arguments print the usage
  and exit with status 2.
- **Chromosome names** must be spelled as the file spells them (`chr1` is not `1`). `info` gives the
  number of chromosomes but not their names; from Java, `chromSizes()` lists them.
- **Memory.** A window is held in memory before it is printed: 8 bytes per base per track. The reader
  also keeps up to 256 MB of the file it has decoded. For a wide window over many tracks, raise the
  heap with `-Xmx`, for example `java -Xmx8g -cp ...`.
- **One window per command.** To read many windows, use the class from Java.

## Using it from Java

```java
import figwig.FigWigReader;

try (FigWigReader fw = new FigWigReader("k562.fw")) {
    System.out.println(fw.genome() + " " + fw.chromSizes().size() + " " + fw.tile());
    System.out.println(fw.tracks().get(0));          // {name=atac.+, dtype=uint8, assay=ATAC-seq, ...}

    String[] chroms = {"chr1", "chr1", "chr2"};
    int[] starts = {1_000_000, 2_500_000, 300_000};

    double[][][] y = fw.read(chroms, starts, 2114);                    // (3 windows, every track, 2114)
    double[][][] atac = fw.read(chroms, starts, 2114,
            fw.select("assay", "ATAC-seq"));                           // the tracks with that attribute
    double[][] one = fw.read(chroms, starts, 2114, 5);                 // track 5 alone: (3, 2114)

    double[] perTile = fw.tileSums(0, "chr1", 0);                      // track 0's sum over each tile
}

// A reader of a few tracks, by name, in the order given.  It acts as if the file held only those.
try (FigWigReader stranded = new FigWigReader("k562.fw", "procap.+", "procap.-")) {
    double[][][] y = stranded.read(new String[]{"chr1"}, new int[]{1_000_000}, 2114);   // (1, 2, 2114)
}
```

- Close the reader when done, as `try` does above. It keeps the file open until then.
- Windows may overlap, repeat and come in any order; each is written to its own row.
- A window that starts before 0, or lies on a chromosome the file does not have, throws
  `IllegalArgumentException`. A file that cannot be read throws `IOException`.
- Track indices passed to `read`, `tileSums` and returned by `select` count the reader's own tracks,
  in the order of `tracks()`.
- One reader may be used from several threads at once.

## Differences from the Python reader

What a user of the Python `FigWigReader` will find different.

### Values

| | Python | Java |
|---|---|---|
| Type of the result | The tracks' own type (`uint8`, `uint16`, `float32`, ...), or one asked for with `dtype=` | Always `double`. Every figwig type fits exactly, so no value changes, but integers are not returned as integers |
| Base past the end of a chromosome | NaN in a float result, 0 in an integer one | Always NaN |
| Chromosome the file does not have | An error, or zeros with `absent="zero"` | Always an error |
| Tile sums | `float32` for tiles, `float64` for zoom levels | `double` for both, holding the same numbers |

### Calls

| Python | Java |
|---|---|
| `fw.read(chroms, starts, width=w, tracks=["a", "b"])`, by name or index | `fw.read(chroms, starts, w, new int[]{i, j})`, by index only. Names are given when the reader is opened |
| `chroms` may be one name for every window | One name per window |
| `fw.select(assay="ATAC-seq")` | `fw.select("assay", "ATAC-seq")`: keys and values in pairs |
| `select` raises when no track has the key | Returns no tracks |
| `select` tells the number `17` from the text `"17"` | Compares as text: both match |
| `FigWigReader(path, tracks={"assay": "ATAC-seq"})` | Not offered. Open with track names, or pass `select(...)` to `read` |
| `fw.stats["sum"]`, a table with named columns | `fw.stats()[i]` is `{sum, nonzero, min, max}` for track `i` |
| `fw.tracks[i]` includes `group`, with attributes in the file's order and of their own types | No `group`; attributes are text and in no particular order |
| `fw.tile_sums()` for every track over the genome, or some tracks and chromosomes | `fw.tileSums(track, chrom, level)`: one track and one chromosome per call |
| `tile_sums` on a track whose sums the file does not keep says so | Throws `IllegalStateException` |
| Holds no open file between reads | Holds the file open until `close()` |

### Not in the Java reader

- Sparse reads (`sparse=True`), sums per label (`labels=`), sums per bin (`bins=`), `sums`, and
  `read_entries`.
- Filling an array the caller supplies (`out=`), choosing the thread count (`n_jobs=`), and reading
  onto a GPU (`device=`).
- The `read_figwig` function, pickling a reader, and the check that the file has not been replaced
  since it was opened.
- Writing, appending, merging, reordering and converting files, and the reader of bigWig files.

### Speed

The two have not been timed against each other. The Python reader decodes many windows and tracks in
one compiled pass. The Java reader reads window by window and track by track, so expect it to be
slower for large batches. It has no compile step on first use; the Python documentation gives about
20 seconds for its first dense read in a new environment, after which the compiled code is cached.
