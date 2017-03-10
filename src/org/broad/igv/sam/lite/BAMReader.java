package org.broad.igv.sam.lite;

import htsjdk.tribble.util.LittleEndianInputStream;
import org.broad.igv.sam.Alignment;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedInputStream;
import java.util.List;

/**
 * Created by jrobinso on 3/9/17.
 */
public class BAMReader {


    int BAM_MAGIC = 21840194;
    int BAI_MAGIC = 21578050;
    char[] SECRET_DECODER = {'=', 'A', 'C', 'x', 'G', 'x', 'x', 'x', 'T', 'x', 'x', 'x', 'x', 'x', 'x', 'N'};
    char[] CIGAR_DECODER = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', '?', '?', '?', '?', '?', '?', '?'};
    byte READ_STRAND_FLAG = 0x10;
    byte MATE_STRAND_FLAG = 0x20;
    int FIRST_OF_PAIR_FLAG = 0x40;
    int SECOND_OF_PAIR_FLAG = 0x80;
    int NOT_PRIMARY_ALIGNMENT_FLAG = 0x100;
    int READ_FAILS_VENDOR_QUALITY_CHECK_FLAG = 0x200;
    int DUPLICATE_READ_FLAG = 0x400;
    int SUPPLEMENTARY_FLAG = 0x800;

    int MAX_GZIP_BLOCK_SIZE = 65536;   //  APPARENTLY.  Where is this documented???
    int DEFAULT_SAMPLING_WINDOW_SIZE = 100;
    int DEFAULT_SAMPLING_DEPTH = 50;
    int MAXIMUM_SAMPLING_DEPTH = 2500;

    BAMIndex bamIndex = null;

    public class BamReader {


    }

   public BAMReader() {

//        this.config = config;
//
//        this.filter = config.filter || new igv.BamFilter();
//
//        this.bamPath = config.url;
//        // Todo - deal with Picard convention.  WHY DOES THERE HAVE TO BE 2?
//        this.baiPath = config.indexURL || this.bamPath + ".bai"; // If there is an indexURL provided, use it!
//        this.headPath = config.headURL || this.bamPath;
//
//
//        this.samplingWindowSize = config.samplingWindowSize == = undefined ? DEFAULT_SAMPLING_WINDOW_SIZE : config.samplingWindowSize;
//        this.samplingDepth = config.samplingDepth == = undefined ? DEFAULT_SAMPLING_DEPTH : config.samplingDepth;
//        if (this.samplingDepth > MAXIMUM_SAMPLING_DEPTH) {
//            igv.log("Warning: attempt to set sampling depth > maximum value of 2500");
//            this.samplingDepth = MAXIMUM_SAMPLING_DEPTH;
//        }
//
//        if (config.viewAsPairs) {
//            this.pairsSupported = true;
//        } else {
//            this.pairsSupported = config.pairsSupported == = undefined ? true : config.pairsSupported;
//        }

    }

    public List<Alignment> readAlignments (String chr, int bpStart, int bpEnd) {

//    function(chr, bpStart, bpEnd) {
//
//        var self = this;
//
//        return new Promise(function(fulfill, reject) {
//
//            getChrIndex(self)
//
//            .
//
//            then(function(chrToIndex) {
//
//                var chrId = chrToIndex[chr],
//
//                        alignmentContainer = new igv.AlignmentContainer(chr, bpStart, bpEnd, self.samplingWindowSize, self.samplingDepth, self.pairsSupported);
//
//                if (chrId == = undefined) {
//                    fulfill(alignmentContainer);
//                } else {
//
//                    getIndex(self)
//                            .then(function(bamIndex) {
//
//                        var chunks = bamIndex.blocksForRange(chrId, bpStart, bpEnd),
//                        promises =[];
//
//
//                        if (!chunks) {
//                            fulfill(null);
//                            reject("Error reading bam index");
//                            return;
//                        }
//                        if (chunks.length == = 0) {
//                            fulfill(alignmentContainer);
//                            return;
//                        }
//
//                        chunks.forEach(function(c) {
//
//                            promises.push(new Promise(function(fulfill, reject) {
//
//                                var fetchMin = c.minv.block,
//                                        fetchMax = c.maxv.block + 65000,   // Make sure we get the whole block.
//                                        range = {start
//                                :fetchMin,size:fetchMax-fetchMin+1
//                            };
//
//                            igvxhr.loadArrayBuffer(self.bamPath,
//                                    {
//                                            headers:self.config.headers,
//                                    range:range,
//                                    withCredentials:self.config.withCredentials
//                            }).then(function(compressed) {
//
//                                var ba = new Uint8Array(igv.unbgzf(compressed)); //new Uint8Array(igv.unbgzf(compressed)); //, c.maxv.block - c.minv.block + 1));
//                                decodeBamRecords(ba, c.minv.offset, alignmentContainer, bpStart, bpEnd, chrId, self.filter);
//
//                                fulfill(alignmentContainer);
//
//                            }).catch(function(obj) {
//                                reject(obj);
//                            });
//
//                            }))
//                        });
//
//
//                        Promise.all(promises).then(function(ignored) {
//                            alignmentContainer.finish();
//                            fulfill(alignmentContainer);
//                        }).catch(function(obj) {
//                            reject(obj);
//                        });
//                    }).catch(reject);
//                }
//            }
//
//            ).catch(reject);

        return null;
        }
//
//
//        function decodeBamRecords (ba, offset, alignmentContainer, min, max, chrId, filter){
//
//            var blockSize,
//                    blockEnd,
//                    alignment,
//                    blocks,
//                    refID,
//                    pos,
//                    bmn,
//                    bin,
//                    mq,
//                    nl,
//                    flag_nc,
//                    flag,
//                    nc,
//                    lseq,
//                    mateRefID,
//                    matePos,
//                    readName,
//                    j,
//                    p,
//                    lengthOnRef,
//                    cigar,
//                    c,
//                    cigarArray,
//                    seq,
//                    seqBytes;
//
//            while (true) {
//
//                blockSize = readInt(ba, offset);
//                blockEnd = offset + blockSize + 4;
//
//                if (blockEnd > ba.length) {
//                    return;
//                }
//
//                alignment = new igv.BamAlignment();
//
//                refID = readInt(ba, offset + 4);
//                pos = readInt(ba, offset + 8);
//
//                if (refID < 0) {
//                    return;   // unmapped reads
//                } else if (refID > chrId || pos > max) {
//                    return;    // off right edge, we're done
//                } else if (refID < chrId) {
//                    continue;   // to left of start, not sure this is possible
//                }
//
//                bmn = readInt(ba, offset + 12);
//                bin = (bmn & 0xffff0000) >> 16;
//                mq = (bmn & 0xff00) >> 8;
//                nl = bmn & 0xff;
//
//                flag_nc = readInt(ba, offset + 16);
//                flag = (flag_nc & 0xffff0000) >> 16;
//                nc = flag_nc & 0xffff;
//
//                alignment.flags = flag;
//                alignment.strand = !(flag & READ_STRAND_FLAG);
//
//                lseq = readInt(ba, offset + 20);
//
//                mateRefID = readInt(ba, offset + 24);
//                matePos = readInt(ba, offset + 28);
//                alignment.fragmentLength = readInt(ba, offset + 32);
//
//                readName = '';
//                for (j = 0; j < nl - 1; ++j) {
//                    readName += String.fromCharCode(ba[offset + 36 + j]);
//                }
//
//                p = offset + 36 + nl;
//
//                lengthOnRef = 0;
//                cigar = '';
//
//
//                cigarArray =[];
//                for (c = 0; c < nc; ++c) {
//                    var cigop = readInt(ba, p);
//                    var opLen = (cigop >> 4);
//                    var opLtr = CIGAR_DECODER[cigop & 0xf];
//                    if (opLtr == 'M' || opLtr == 'EQ' || opLtr == 'X' || opLtr == 'D' || opLtr == 'N' || opLtr == '=')
//                        lengthOnRef += opLen;
//                    cigar = cigar + opLen + opLtr;
//                    p += 4;
//
//                    cigarArray.push({len:opLen, ltr:opLtr});
//                }
//                alignment.cigar = cigar;
//                alignment.lengthOnRef = lengthOnRef;
//
//                if (alignment.start + alignment.lengthOnRef < min)
//                    continue;  // Record out-of-range "to the left", skip to next one
//
//
//                seq = '';
//                seqBytes = (lseq + 1) >> 1;
//                for (j = 0; j < seqBytes; ++j) {
//                    var sb = ba[p + j];
//                    seq += SECRET_DECODER[(sb & 0xf0) >> 4];
//                    seq += SECRET_DECODER[(sb & 0x0f)];
//                }
//                seq = seq.substring(0, lseq);  // seq might have one extra character (if lseq is an odd number)
//
//                p += seqBytes;
//                alignment.seq = seq;
//
//
//                if (lseq == = 1 && String.fromCharCode(ba[p + j] + 33) == = "*") {
//                    // TODO == how to represent this?
//                } else {
//                    alignment.qual =[];
//                    for (j = 0; j < lseq; ++j) {
//                        alignment.qual.push(ba[p + j]);
//                    }
//                }
//                p += lseq;
//
//
//                alignment.start = pos;
//                alignment.mq = mq;
//                alignment.readName = readName;
//                alignment.chr = self.indexToChr[refID];
//
//                if (mateRefID >= 0) {
//                    alignment.mate = {
//                            chr:self.indexToChr[mateRefID],
//                            position:matePos,
//                            strand:!(flag & MATE_STRAND_FLAG)
//                    };
//                }
//
//
//                alignment.tagBA = new Uint8Array(ba.buffer.slice(p, blockEnd));  // decode thiese on demand
//                p += blockEnd;
//
//                if (!min || alignment.start <= max &&
//                        alignment.start + alignment.lengthOnRef >= min &&
//                        filter.pass(alignment)) {
//                    if (chrId == = undefined || refID == chrId) {
//                        blocks = makeBlocks(alignment, cigarArray);
//                        alignment.blocks = blocks.blocks;
//                        alignment.insertions = blocks.insertions;
//                        alignmentContainer.push(alignment);
//                    }
//                }
//                offset = blockEnd;
//            }
//            // Exits via top of loop.
//        }
//
//        /**
//         * Split the alignment record into blocks as specified in the cigarArray.  Each aligned block contains
//         * its portion of the read sequence and base quality strings.  A read sequence or base quality string
//         * of "*" indicates the value is not recorded.  In all other cases the length of the block sequence (block.seq)
//         * and quality string (block.qual) must == the block length.
//         *
//         * NOTE: Insertions are not yet treated // TODO
//         *
//         * @param record
//         * @param cigarArray
//         * @returns array of blocks
//         */
//        function makeBlocks (record, cigarArray){
//
//            var blocks =[],
//            insertions,
//                    seqOffset = 0,
//                    pos = record.start,
//                    len = cigarArray.length,
//                    blockSeq,
//                    blockQuals,
//                    gapType,
//                    minQ = 5,  //prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MIN)
//                    maxQ = 20; //prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MAX)
//
//            for (var i = 0; i < len; i++) {
//
//                var c = cigarArray[i];
//
//                switch (c.ltr) {
//                    case 'H':
//                        break; // ignore hard clips
//                    case 'P':
//                        break; // ignore pads
//                    case 'S':
//                        seqOffset += c.len;
//                        gapType = 'S';
//                        break; // soft clip read bases
//                    case 'N':
//                        pos += c.len;
//                        gapType = 'N';
//                        break;  // reference skip
//                    case 'D':
//                        pos += c.len;
//                        gapType = 'D';
//                        break;
//                    case 'I':
//                        blockSeq = record.seq == = "*" ? "*" : record.seq.substr(seqOffset, c.len);
//                        blockQuals = record.qual ? record.qual.slice(seqOffset, c.len) : undefined;
//                        if (insertions == = undefined) insertions =[];
//                        insertions.push({start:pos, len:c.len, seq:blockSeq, qual:blockQuals});
//                        seqOffset += c.len;
//                        break;
//                    case 'M':
//                    case 'EQ':
//                    case '=':
//                    case 'X':
//
//                        blockSeq = record.seq == = "*" ? "*" : record.seq.substr(seqOffset, c.len);
//                        blockQuals = record.qual ? record.qual.slice(seqOffset, c.len) : undefined;
//                        blocks.push({start:pos, len:c.len, seq:blockSeq, qual:blockQuals, gapType:gapType});
//                        seqOffset += c.len;
//                        pos += c.len;
//
//                        break;
//
//                    default:
//                        console.log("Error processing cigar element: " + c.len + c.ltr);
//                }
//            }
//
//            return {blocks:blocks, insertions:insertions};
//
//        }
//    }


//    void readHeader() {
//
//                long len = bamIndex.firstAlignmentBlock + MAX_GZIP_BLOCK_SIZE;   // Insure we get the complete compressed block containing the header
//
//        //LittleEndianInputStream parser = new LittleEndianInputStream(new BufferedInputStream(ParsingUtils.openInputStream(indexURL)));

//
//        igvxhr.loadArrayBuffer(self.bamPath,
//                        {
//                                headers:self.config.headers,
//
//                        range:{
//                    start:
//                    0, size:len
//                },
//
//                withCredentials:
//                self.config.withCredentials
//                }).then(function(compressedBuffer) {
//
//                    var unc = igv.unbgzf(compressedBuffer, len),
//                            uncba = new Uint8Array(unc),
//                            magic = readInt(uncba, 0),
//                            samHeaderLen = readInt(uncba, 4),
//                            samHeader = '',
//                            genome = igv.browser ? igv.browser.genome : null;
//
//                    for (var i = 0; i < samHeaderLen; ++i) {
//                        samHeader += String.fromCharCode(uncba[i + 8]);
//                    }
//
//                    var nRef = readInt(uncba, samHeaderLen + 8);
//                    var p = samHeaderLen + 12;
//
//                    self.chrToIndex = {};
//                    self.indexToChr =[];
//                    for (var i = 0; i < nRef; ++i) {
//                        var lName = readInt(uncba, p);
//                        var name = '';
//                        for (var j = 0; j < lName - 1; ++j) {
//                            name += String.fromCharCode(uncba[p + 4 + j]);
//                        }
//                        var lRef = readInt(uncba, p + lName + 4);
//                        //dlog(name + ': ' + lRef);
//
//                        if (genome && genome.getChromosomeName) {
//                            name = genome.getChromosomeName(name);
//                        }
//
//                        self.chrToIndex[name] = i;
//                        self.indexToChr.push(name);
//
//                        p = p + 8 + lName;
//                    }
//
//                    fulfill();
//
//                }).catch(reject);
//
//
//    }

//    //
//    function getIndex(bam) {
//
//        return new Promise(function(fulfill, reject) {
//
//            if(bam.index)
//
//            {
//                fulfill(bam.index);
//            }
//
//            else
//
//            {
//                igv
//                        .loadBamIndex(bam.baiPath, bam.config)
//                        .then(function(index) {
//                bam.index = index;
//                fulfill(bam.index);
//            })
//                .catch(reject);
//            }
//        });
//    }
//
//
//    function getChrIndex(bam) {
//
//        return new Promise(function(fulfill, reject) {
//
//            if(bam.chrToIndex)
//
//            {
//                fulfill(bam.chrToIndex);
//            }
//
//            else
//
//            {
//                bam.readHeader().then(function() {
//                fulfill(bam.chrToIndex);
//            }).catch(reject);
//            }
//        });
//    }

//    function readInt(ba, offset) {
//        return (ba[offset + 3] << 24) | (ba[offset + 2] << 16) | (ba[offset + 1] << 8) | (ba[offset]);
//    }
//
//    function readShort(ba, offset) {
//        return (ba[offset + 1] << 8) | (ba[offset]);
//    }


}
