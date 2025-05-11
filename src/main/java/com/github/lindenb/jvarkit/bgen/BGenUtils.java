/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.bgen;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.math.BigInteger;
import java.nio.charset.Charset;
import java.util.BitSet;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.math.MathUtils;

import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.tribble.TribbleException;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class BGenUtils {
	public static final VCFFormatHeaderLine GP_format_header_line = new VCFFormatHeaderLine(
			"GP",VCFHeaderLineCount.G,VCFHeaderLineType.Float,"Genotype call probabilities");
	public static final VCFFormatHeaderLine HP_format_header_line = new VCFFormatHeaderLine(
			"HP",VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.Float,"Haplotype call probabilities");
	public static final VCFInfoHeaderLine OFFSET_format_header_line = new VCFInfoHeaderLine(
			"OFFSET",1,VCFHeaderLineType.String,"Physical offset in BGEN file; Encoded as String because it can be a long/int64.");

	
	public enum Compression {
		e_NoCompression,
		e_ZlibCompression,
		e_ZstdCompression
		};
		
	public enum Layout {
		e_Layout1,
		e_Layout2
		};
	
	BGenUtils() {
		
	}
	
	static final Charset ENCODING=Charset.forName("UTF-8");
	static final byte[] BGEN_MAGIC = "bgen".getBytes(ENCODING);
	static final String FILE_SUFFIX = ".bgen";
	static final String VCF_INDEX_SUFFIX = ".bgenix"+FileExtensions.COMPRESSED_VCF;

	/** write 'N' bytes */
	static void writeNBytes(BinaryCodec binaryCodec,long n,byte value) {
		while(n>0) {
			binaryCodec.writeByte(value);
			n--;
			}
		}
	
	/** write 'N' null bytes */
	static void writeNBytes(BinaryCodec binaryCodec,long n) {
		writeNBytes(binaryCodec,n,(byte)0);
		}

	
    /** convert a long to int, throws a tribble exception if the number is greater than Integer.MAX_VALUE */
    static int longToUnsignedInt(final long n) {
        if(n>Integer.MAX_VALUE || n<0L) throw new TribbleException("Cannot convert "+n+"L to int");
        return (int)n;
    	}
    /** read the length of a string as unsigned short and read the string */
    static String readStringUInt16(final BinaryCodec binaryCodec) throws IOException {
    	return readString(binaryCodec, binaryCodec.readUShort());
    	}
    
    static void writeStringUInt16(final BinaryCodec codec,final String s) throws IOException {
    	final byte[] b=s.getBytes(ENCODING);
    	codec.writeUShort(b.length);
    	codec.writeBytes(b);
    	}
    
    static void writeStringUInt32(final BinaryCodec codec,final String s) throws IOException {
    	final byte[] b=s.getBytes(ENCODING);
    	codec.writeUInt(b.length);
    	codec.writeBytes(b);
    	}
    
    /** read the length of a string as unsigned short and read the string */
    static String readStringUInt32(final BinaryCodec binaryCodec) throws IOException {
    	return readString(binaryCodec, longToUnsignedInt(binaryCodec.readUInt()));
    	}
    
    /** problem with bitSet, when using toByteArray, it goes to the last SET bit, this function use padding */
    static byte[] toByteArray(final BitSet bitSet,final int n) {
    	// non, pas sur du size()  . if(bitSet.size()>8*n) throw new IllegalStateException("the size of the bitset="+bitSet.size()+" is larger than the number of 8*bytes="+n);
    	final byte[] raw = bitSet.toByteArray();
    	if(raw.length>n) throw new IllegalStateException();
		if(raw.length==n) return raw;
		final byte[] padding = new byte[n];
		System.arraycopy(raw, 0, padding, 0,raw.length);
		return padding;
    }
    
    /**
     * read a string of length 'len'
     */
   private static String readString(final BinaryCodec binaryCodec, final int len) throws IOException {
        return new String(readNBytes(binaryCodec, len),ENCODING);
    	}
    
    /** fully read 'len' bytes */
   static byte[] readNBytes(final BinaryCodec binaryCodec, final int len) throws IOException {
    	if(len<0) throw new IllegalArgumentException("negative bytes length:" + len);
        final byte[] bytes = new byte[len];
        binaryCodec.readBytes(bytes);
        return bytes;
    	}

   
   static int calculateTotalCombinationsForLayout2(boolean phased,final int n_alleles,final int m_ploidy) {
		final BigInteger bi;
		if(phased) {
			bi = BigInteger.valueOf(n_alleles).pow(m_ploidy);
			}
		else
			{
			int n = n_alleles;
		    int k = m_ploidy;
			bi = MathUtils.factorial(n + k - 1).divide(MathUtils.factorial(k).multiply(MathUtils.factorial(n - 1)));
			}
		return bi.intValueExact();
		}
   
   	/** ByteArrayOutputStream but we expose the internal buffer */
    protected static class ByteBuffer extends ByteArrayOutputStream {
		byte[] getBuffer() { return super.buf;}
		public void copyTo(InputStream in,int nbytes) throws IOException {
			IOUtils.copyTo(in,this,nbytes);
			}
		ByteArrayInputStream toByteArrayInputStream() {
			return new ByteArrayInputStream(getBuffer(), 0, size());
			}
		}
    
    
    /** bit reader */
	static class BitReader  {
		private final InputStream in;
		private byte curr;
		private int offset = 8;
		private final long max_to_read;
		private long nRead = 0L;
		BitReader(final InputStream in,long max_to_read) {
			this.in = in;
			this.max_to_read=max_to_read;
			}
		BitReader(final InputStream in) {
			this(in,-1L);
			}
		/** return 0/1 or -1 */
		int read() throws IOException {
			if(max_to_read!=-1L && nRead>=max_to_read) return -1;
			if (offset >= 8) {
				int c = in.read();
				if (c == -1)
					return -1;
				this.curr = (byte) c;
				offset = 0;
				}
			final int bit = (curr >> (offset)) & 1;
			offset++;
			nRead++;
			return bit;
			}
	
			
		@Override
		public String toString() {
			return "BitReader("+(int)curr+","+offset+")";
			}
		}
	
	
	
    /** bit-based number reader reader */
	static class BitNumReader  {
		private final double max;
		private final BitReader delegate;
		private final int nbits;
		BitNumReader(BitReader delegate,int nbits) {
			if (nbits <= 0 || nbits > 32) {
			    throw new IllegalArgumentException("The number of bits (nbits) must be between 1 and 32.");
				}
			this.max = Math.pow(2,nbits)-1;
			this.delegate = delegate;
			this.nbits = nbits;
			}
		public double next() throws IOException {
			int t=0;
			int n = this.nbits;
			while(n>0) {
				final int c = delegate.read();
				if(c==-1) throw new IOException("cannot read "+n+"th bit");
				t = (t << 1) | c;
				// t=t*2 + c;
				--n;
				}
			return t/max;
			}
		}
	

	static class BitWriter {
	    private final OutputStream out;
	    private byte curr = '\0'; // Current byte being written
	    private int offset = 0; // Tracks the position of the next bit in the byte

	    BitWriter(final OutputStream out) {
	        this.out = out;
	    	}

	    void write(boolean bit) throws IOException {
	        // Set or clear the bit at the current offset
	        if (bit) {
	            curr |= (1 << offset); // Set the bit
	        	}
	        offset++;

	        // If the byte is full (8 bits written), write it to the OutputStream
	        if (offset == 8) flushByte();
	    	}
	    
	    private void flushByte() throws IOException {
	    	out.write(curr);
            curr = '\0'; // Reset the current byte buffer
            offset = 0; // Reset the offset
	    	}
	    
	    void flush() throws IOException {
	        // If there are remaining bits, flush them as a partial byte
	        if (offset > 0) flushByte();
	        out.flush(); // Ensure the OutputStream is flushed
	    }

	}
	
	public static class BitNumWriter {
	    private final BitWriter delegate;
	    private final int nbits;
	    private final double max;

	    /**
	     * Constructor for BitNumWriter
	     * @param delegate The BitWriter instance used for writing bits
	     * @param nbits The number of bits to represent each number
	     */
	    public BitNumWriter(BitWriter delegate, int nbits) {
	        if (nbits <= 0 || nbits > 32) {
	            throw new IllegalArgumentException("The number of bits (nbits) must be between 1 and 32.");
	        }
	        this.delegate = delegate;
	        this.nbits = nbits;
	        this.max = Math.pow(2, nbits) - 1; // Maximum value representable with nbits
	    }

	    /**
	     * Writes a normalized number (between 0 and 1) as nbits.
	     * @param value The number to write (must be in the range [0, 1])
	     * @throws IOException If writing to the BitWriter fails
	     */
	    public void write(double value) throws IOException {
	        if (value < 0.0 || value > 1.0) {
	            throw new IllegalArgumentException("The value "+value+" must be in the range [0, 1].");
	        }

	        // Scale the value to an integer representation
	        final int scaledValue = (int) Math.round(value * max);

	        // Write the scaled value bit by bit
	        for (int i = nbits - 1; i >= 0; i--) {
	            boolean bit = (scaledValue & (1 << i)) != 0;
	            delegate.write(bit);
	        }
	    }


	}
}
