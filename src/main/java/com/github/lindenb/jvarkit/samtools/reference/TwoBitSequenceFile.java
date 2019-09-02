/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
package com.github.lindenb.jvarkit.samtools.reference;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;

public class TwoBitSequenceFile implements ReferenceSequenceFile {
    /** standard suffix of 2bit files */
    public static final String SUFFIX = ".2bit";
    // https://github.com/rpique/UCSC-Browser-code-add-ons/blob/ba859c047e46d8074b93700d34cb601fa9ba4288/src/lib/dnautil.c
    private static final int MASKED_BASE_BIT = 8;
    // Numerical values for bases.
    private static final int T_BASE_VAL = 0;
    // private static final int U_BASE_VAL =0;
    private static final int C_BASE_VAL = 1;
    private static final int A_BASE_VAL = 2;
    private static final int G_BASE_VAL = 3;
    private static final int N_BASE_VAL = 4;// Used in 1/2 byte representation.

    private static final byte valToNucl[] = new byte[(N_BASE_VAL | MASKED_BASE_BIT) + 1];

    {
        {
            {
                valToNucl[T_BASE_VAL] = valToNucl[T_BASE_VAL | MASKED_BASE_BIT] = 't';
                valToNucl[C_BASE_VAL] = valToNucl[C_BASE_VAL | MASKED_BASE_BIT] = 'c';
                valToNucl[A_BASE_VAL] = valToNucl[A_BASE_VAL | MASKED_BASE_BIT] = 'a';
                valToNucl[G_BASE_VAL] = valToNucl[G_BASE_VAL | MASKED_BASE_BIT] = 'g';
                valToNucl[N_BASE_VAL] = valToNucl[N_BASE_VAL | MASKED_BASE_BIT] = 'n';
            }
        }
    }

    /*
     * Signature into 2bit file (2 bits per nucleotide DNA file) plus information on N and masked
     * bases.
     */
    private final static int twoBitSig = 0x1A412743;

    /* Signature of byte-swapped two-bit file. */
    private final static int twoBitSwapSig = 0x4327411A;

    private final SeekableStream seekableStream;
    private final ByteOrder byteOrder;
    private Map<String, TwoBitIndex> seq2index;
    private final SAMSequenceDictionary dict;
    /** iterator for {@link #nextSequence()} */
    private Iterator<SAMSequenceRecord> entryIterator = null;
    /*
     * Cache information about last sequence accessed, including nBlock and mask block. This doesn't
     * include the data. This speeds fragment reads.
     */
    private TwoBit seqCache = null;

    private static class Block {
        int count;
        int starts[];
        int sizes[];
    }

    private static class TwoBit {
        String name; /* Name of sequence. */
        @SuppressWarnings("unused")
        byte[] data; /* DNA at two bits per base. */
        int size; /* Size of this sequence. */
        Block nBlock;
        // int nBlockCount; /* Count of blocks of Ns. */
        // int[] nStarts; /* Starts of blocks of Ns. */
        // int[] nSizes; /* Sizes of blocks of Ns. */
        Block maskBlock;
        // int maskBlockCount; /* Count of masked blocks. */
        // int[] maskStarts; /* Starts of masked regions. */
        // int[] maskSizes; /* Sizes of masked regions. */
        @SuppressWarnings("unused")
        int reserved; /* Reserved for future expansion. */

        /* not in original C-Struct file. Offset of data for seqCache sequence */
        long dataOffsetCache = -1L;
    }

    private static class TwoBitIndex {
        String name;
        int seqIndex;
        long offset;
    }

    public TwoBitSequenceFile(final Path path) throws IOException {
        this(path, true);
    }

    public TwoBitSequenceFile(final Path path, final boolean truncateNamesAtWhitespace)
            throws IOException {
    		this(path.toString(),truncateNamesAtWhitespace);
    		}
    
    public TwoBitSequenceFile(final String path, final boolean truncateNamesAtWhitespace)
            throws IOException {
    		this(SeekableStreamFactory.getInstance().getStreamFor(path),truncateNamesAtWhitespace);
    		}
    
    public TwoBitSequenceFile(final SeekableStream seekableStream, final boolean truncateNamesAtWhitespace) throws IOException {
    	this.seekableStream = seekableStream;
    	
    	final byte array[]= new byte[Integer.BYTES];
    	this.seekableStream.readFully(array);
        final int sig =  ByteBuffer.wrap(array).getInt();
        if (sig == twoBitSig) {
            this.byteOrder = ByteOrder.BIG_ENDIAN;
        } else if (sig == twoBitSwapSig) {
            this.byteOrder = ByteOrder.LITTLE_ENDIAN;
        } else {
            throw new IOException("Cannot read header from " + seekableStream.getSource());
        }
        final int version = this.readInt();
        if (version != 0) {
            throw new IOException(
                    "Can only handle version 0 or version 1 of this file. This is version "
                            + version);
        }
        final int seqCount = this.readInt();
        /* int reserved ignored */ this.readInt();
        this.seq2index = new HashMap<>(seqCount);
        for (int i = 0; i < seqCount; i++) {
            String seqName = this.readString();
            if (truncateNamesAtWhitespace) {
                int ws = 0;
                for (ws = 0; ws < seqName.length(); ++ws) {
                    if (Character.isWhitespace(seqName.charAt(ws)))
                        break;
                }
                seqName = seqName.substring(0, ws);
            }
            if (this.seq2index.containsKey(seqName)) {
                throw new IOException("duplicate sequence name \"" + seqName + "\" in " + seekableStream.getSource());
            }
            final TwoBitIndex twoBitIndex = new TwoBitIndex();
            twoBitIndex.name = seqName;
            twoBitIndex.seqIndex = this.seq2index.size();
            twoBitIndex.offset = this.readInt();
            this.seq2index.put(twoBitIndex.name, twoBitIndex);
        }

        final List<SAMSequenceRecord> ssrs = new ArrayList<>(this.seq2index.size());

        for (final Iterator<String> iter = this.seq2index.values().stream()
                .sorted((A, B) -> Integer.compare(A.seqIndex, B.seqIndex)).map(R -> R.name)
                .iterator(); iter.hasNext();) {
            final String contig = iter.next();
            final TwoBit tbi = getTwoBitSeqHeader(contig);
            final int length = tbi.size;
            ssrs.add(new SAMSequenceRecord(contig, length));
        }
        this.dict = new SAMSequenceDictionary(ssrs);
        this.entryIterator = this.dict.getSequences().iterator();

    }

    private byte[] query(final Locatable loc, boolean doMask) throws IOException {
        int remainder, midStart, midEnd;
        final TwoBitIndex tbi = this.seq2index.get(loc.getContig());
        if (tbi == null) {
            throw new IllegalArgumentException("cannot find sequence " + loc.getContig());
        }
        final TwoBit twoBit = getTwoBitSeqHeader(loc.getContig());
        int fragStart = loc.getStart() - 1;
        int fragEnd = loc.getEnd();

        /* validate range. */
        if (fragEnd > twoBit.size) {
            throw new IllegalArgumentException("twoBitReadSeqFrag in " + loc.getContig() + " end ("
                    + fragEnd + ") >= seqSize (" + twoBit.size + ")");
        }
        final int outSize = fragEnd - fragStart;
        if (outSize < 1) {
            throw new IllegalArgumentException("twoBitReadSeqFrag in  " + loc.getContig()
                    + " start (" + fragStart + ") >= end (" + fragEnd + ")");
        }

        int packedStart = (fragStart >> 2);
        int packedEnd = ((fragEnd + 3) >> 2);
        int packByteCount = packedEnd - packedStart;
        this.seekableStream.seek(this.seekableStream.position() + packedStart);
        final byte packed[] = new byte[packByteCount];
        this.seekableStream.readFully(packed);
        final byte dna[] = new byte[outSize];
        int dna_idx = 0;
        int packed_idx = 0;

        /* Handle case where everything is in one packed byte */
        if (packByteCount == 1) {
            int pOff = (packedStart << 2);
            int pStart = fragStart - pOff;
            int pEnd = fragEnd - pOff;
            int partial = Byte.toUnsignedInt(packed[0]);
            assert(pEnd <= 4);
            assert(pStart >= 0);
            for (int i = pStart; i < pEnd; ++i) {
                dna[dna_idx++] = valToNt((partial >> (6 - i - i)) & 3);
            }
        } else {
            /* Handle partial first packed byte. */
            midStart = fragStart;
            remainder = (fragStart & 3);
            if (remainder > 0) {
                int partial = Byte.toUnsignedInt(packed[packed_idx++]);
                int partCount = 4 - remainder;
                for (int i = partCount - 1; i >= 0; --i) {
                    dna[dna_idx + i] = valToNt(partial & 3);
                    partial >>= 2;
                }
                midStart += partCount;
                dna_idx += partCount;
            }

            /* Handle middle bytes. */
            remainder = fragEnd & 3;
            midEnd = fragEnd - remainder;

            for (int i = midStart; i < midEnd; i += 4) {
                int b = Byte.toUnsignedInt(packed[packed_idx++]);
                dna[dna_idx + 3] = valToNt(b & 3);
                b >>= 2;
                dna[dna_idx + 2] = valToNt(b & 3);
                b >>= 2;
                dna[dna_idx + 1] = valToNt(b & 3);
                b >>= 2;
                dna[dna_idx + 0] = valToNt(b & 3);
                dna_idx += 4;
            }

            if (remainder > 0) {
                int part = Byte.toUnsignedInt(packed[packed_idx]);
                part >>= (8 - remainder - remainder);
                for (int i = remainder - 1; i >= 0; --i) {
                    dna[dna_idx + i] = valToNt(part & 3);
                    part >>= 2;
                }
            }
        }

        for (int side = 0; side < 2; ++side) {
            final Block block = (side == 0 ? twoBit.nBlock : twoBit.maskBlock);
            if (block.count == 0)
                continue;
            if (side == 1) {
                if (!doMask)
                    continue;
                toUpperN(dna);
            }

            int startIx = findGreatestLowerBound(block.count, block.starts, fragStart);
            for (int i = startIx; i < block.count; ++i) {
                int s = block.starts[i];
                int e = s + block.sizes[i];
                if (s >= fragEnd)
                    break;
                if (s < fragStart)
                    s = fragStart;
                if (e > fragEnd)
                    e = fragEnd;
                if (s < e) {
                    final int arrayStart = s - fragStart;
                    final int arrayLen = e - s;
                    if (side == 0) {
                        Arrays.fill(dna, arrayStart, arrayStart + arrayLen, (byte) 'n');
                        // memset(seq->dna + s - fragStart, 'n', e - s);
                    } else {
                        for (int x = 0; x < arrayLen; ++x) {
                            dna[arrayStart + x] =  (byte) Character.toLowerCase(dna[arrayStart + x]);
                        }
                    }
                }
            }
        }
        return dna;
    }

    @Override
    public SAMSequenceDictionary getSequenceDictionary() {
        return this.dict;
    }

    @Override
    public ReferenceSequence getSubsequenceAt(final String contig, long start, long stop) {
        if (start > Integer.MAX_VALUE)
            throw new SAMException("start is too large: " + stop);
        if (stop > Integer.MAX_VALUE)
            throw new SAMException("stop is too large: " + stop);
        if (this.dict.getSequence(contig) == null)
            return null;
        try {
            final byte bases[] = query(new Interval(contig, (int) start, (int) stop), false);
            return new ReferenceSequence(contig, (int) start - 1 /* 0-based */, bases);
        } catch (IOException err) {
            throw new RuntimeIOException(err);
        }
    }

    @Override
    public ReferenceSequence getSequence(final String contig) {
        final SAMSequenceRecord ssr = this.dict.getSequence(contig);
        if (ssr == null)
            return null;
        return getSubsequenceAt(contig, 1, ssr.getSequenceLength());
    }

    @Override
    public boolean isIndexed() {
        return true;
    }

    @Override
    public ReferenceSequence nextSequence() {
        if (!this.entryIterator.hasNext()) {
            return null;
        }
        return getSequence(this.entryIterator.next().getSequenceName());
    }

    @Override
    public void reset() {
        this.entryIterator = this.dict.getSequences().iterator();
    }

    private byte valToNt(int b) {
        /*
         * #define MASKED_BASE_BIT 8
         * 
         * Numerical values for bases. #define T_BASE_VAL 0 #define U_BASE_VAL 0 #define C_BASE_VAL
         * 1 #define A_BASE_VAL 2 #define G_BASE_VAL 3 #define N_BASE_VAL 4 Used in 1/2 byte
         * representation.
         * 
         * valToNt[T_BASE_VAL] = valToNt[T_BASE_VAL|MASKED_BASE_BIT] = 't'; valToNt[C_BASE_VAL] =
         * valToNt[C_BASE_VAL|MASKED_BASE_BIT] = 'c'; valToNt[A_BASE_VAL] =
         * valToNt[A_BASE_VAL|MASKED_BASE_BIT] = 'a'; valToNt[G_BASE_VAL] =
         * valToNt[G_BASE_VAL|MASKED_BASE_BIT] = 'g'; valToNt[N_BASE_VAL] =
         * valToNt[N_BASE_VAL|MASKED_BASE_BIT] = 'n';
         */
        return valToNucl[b];
    }

    private ByteBuffer mustRead(int nBytes) throws IOException {
    	final byte array[] = new byte[nBytes];
    	this.seekableStream.readFully(array);
        return ByteBuffer.wrap(array).order(this.byteOrder);
    }

    private int readInt() throws IOException {
        return mustRead(Integer.BYTES).getInt();
    }

    private byte readByte() throws IOException {
        return mustRead(Byte.BYTES).get();
    }

    private String readString() throws IOException {
        final int nchar = Byte.toUnsignedInt(this.readByte());
        final byte array[] = new byte[nchar*Byte.BYTES];
        this.seekableStream.readFully(array);
        return new String(array);
    }

    @Override
    public void close() throws IOException {
        this.seekableStream.close();
    }


    private TwoBit getTwoBitSeqHeader(final String name) throws IOException
    /*
     * get the sequence header information using the cache. Position file right at data.
     */
    {
        if (this.seqCache != null && this.seqCache.name.equals(name)) {
            this.seekableStream.seek(this.seqCache.dataOffsetCache);
        } else {
            // fetch new and cache
            this.seqCache = readTwoBitSeqHeader(name);
        }
        return this.seqCache;
    }



    private TwoBit readTwoBitSeqHeader(final String name) throws IOException
    /*
     * read a sequence header, nBlocks and maskBlocks from a twoBit file, leaving file pointer at
     * data block
     */
    {
        final TwoBit twoBit = new TwoBit();
        twoBit.name = name;


        /* Find offset in index and seek to it */
        twoBitSeekTo(name);

        /* Read in seqSize. */
        twoBit.size = readInt();

        /* Read in blocks of N. */
        twoBit.nBlock = readBlockCoords();


        /* Read in masked blocks. */
        twoBit.maskBlock = readBlockCoords();

        /* Reserved word. */
        twoBit.reserved = readInt();

        twoBit.dataOffsetCache = this.seekableStream.position();

        return twoBit;
    }

    private void toUpperN(byte array[]) {

    }



    private Block readBlockCoords() throws IOException
    /*
     * Read in blockCount, starts and sizes from file. (Same structure used for both blocks of N's
     * and masked blocks.)
     */
    {
        final Block block = new Block();
        block.count = readInt();
        if (block.count == 0) {
            block.starts = null;
            block.sizes = null;
        } else {
            block.starts = new int[block.count];
            block.sizes = new int[block.count];
            final byte array[] = new byte[block.count * Integer.BYTES];
            
            ByteBuffer buffer =ByteBuffer.wrap(array).order(this.byteOrder);
            this.seekableStream.readFully(array);
            block.starts = new int[block.count];
            for(int i=0;i<block.count;++i )
            	{
            	block.starts[i]=buffer.getInt();
            	}
            
            buffer =ByteBuffer.wrap(array).order(this.byteOrder);
            this.seekableStream.readFully(array);
            block.sizes = new int[block.count];
            for(int i=0;i<block.count;++i )
            	{
            	block.sizes[i]=buffer.getInt();
            	}
        }
        return block;
    }


    /* Seek to start of named record. Abort if can't find it. */
    private void twoBitSeekTo(final String name) throws IOException {
        final TwoBitIndex index = this.seq2index.get(name);
        if (index == null)
            throw new IllegalArgumentException("not int dict" + name);
        this.seekableStream.seek(index.offset);
    }

    private int findGreatestLowerBound(int blockCount, int[] pos, int val)
    /*
     * Find index of greatest element in posArray that is less than or equal to val using a binary
     * search.
     */
    {
        int startIx = 0, endIx = blockCount - 1, midIx;
        int posVal;

        for (;;) {
            if (startIx == endIx) {
                posVal = pos[startIx];
                if (posVal <= val || startIx == 0)
                    return startIx;
                else
                    return startIx - 1;
            }
            midIx = ((startIx + endIx) >> 1);
            posVal = pos[midIx];
            if (posVal < val)
                startIx = midIx + 1;
            else
                endIx = midIx;
        }
    }

    @Override
    public String toString() {
        return "TwoBitSequenceFile(" +this.seekableStream.getSource() + ")";
    }

    public static void main(final String[] args) {
        try {
            TwoBitSequenceFile r = new TwoBitSequenceFile(Paths
                    .get("/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.2bit"));
            byte array[] = r.query(new Interval("RF01", 1, 20), false);
            System.err.println(new String(array));
            System.err.println(r.getSequence("RF11").getBaseString());
            r.close();
            System.err.println("done");
        } catch (Exception err) {
            err.printStackTrace();
        }
    }
}
