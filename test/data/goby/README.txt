These assist with TestIteratedSortedAlignment2.

The index is "small-synth.fa". The reads are "seq-var-reads.fa".

To create the gsnap index, first create compact-reads from the files

   java -jar goby.jar -m fasta-to-compact -d -o small-synth.compact-reads small-synth.fa
   java -jar goby.jar -m fasta-to-compact -d -o seq-var-reads.compact-reads seq-var-reads.fa

Create the GSNAP index

   mkdir SMALL_SYNTH_DB-gsnap
   cd SMALL_SYNTH_DB-gsnap
   ../../gmap-icb/util/gmap_setup -d index -B ../../gmap-icb/util ../small-synth.fa
   make -f Makefile.index coords
   make -f Makefile.index gmapdb
   # The following augments the index for Bisulfite, not necessary
   ../../gmap-icb/src/cmetindex -d index -D . -F .
   cd ..

Align the output. To make the output useful for debugging, make sure to
compile goby-ccp with debugging enabled by editing the top of C_Alignments.cc
and set
    #define C_WRITE_API_WRITE_ALIGNMENT_DEBUG

Align command:

   # unthreaded with 6 mismatches allowed
   ../gmap-icb/src/gsnap -t 0 -m 6 -A goby --goby-output seq-var-reads-gsnap -D SMALL_SYNTH_DB-gsnap -d index seq-var-reads.compact-reads > seq-var-reads-gsnap.exec-output.txt

Finally, since the iterator we are testing requires a SORTED alignment run

   java -jar goby.jar -m sort -d -o sorted-seq-var-reads-gsnap seq-var-reads-gsnap.entries

The files named *-hybrid have been compressed with the Hybrid-1 codec of Goby 2.0 with the following options:

goby 1g ca GDFQPGI-pickrellNA18486_yale -o GDFQPGI-pickrellNA18486_yale-hybrid -x MessageChunksWriter:codec=hybrid-1 -x MessageChunksWriter:template-compression=true -x AlignmentCollectionHandler:enable-domain-optimizations=true -x AlignmentWriterImpl:permutate-query-indices=true -x AlignmentCollectionHandler:ignore-read-origin=true -x AlignmentCollectionHandler:debug-level=1

