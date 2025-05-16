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
import java.io.FilterInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigInteger;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.function.Consumer;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.math.MathUtils;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.seekablestream.SeekableBufferedStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.LocationAware;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class BGenUtils {
	private static final Logger LOG = Logger.of(BGenUtils.class).setDebug();

	public static final VCFFormatHeaderLine GP_format_header_line = new VCFFormatHeaderLine(
			"GP",VCFHeaderLineCount.G,VCFHeaderLineType.Float,"Genotype call probabilities");
	public static final VCFFormatHeaderLine HP_format_header_line = new VCFFormatHeaderLine(
			"HP",VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.Float,"Haplotype call probabilities");
	public static final VCFInfoHeaderLine OFFSET_format_header_line = new VCFInfoHeaderLine(
			"OFFSET",1,VCFHeaderLineType.String,"Physical offset in BGEN file; Encoded as String because it can be a long/int64.");

	
	public enum Compression {
		NONE,
		ZLIB,
		ZSTD
		};
		
	public enum Layout {
		LAYOUT_1,
		LAYOUT_2
		};
	
	BGenUtils() {
		
	}
	
	public static final Charset ENCODING=Charset.forName("UTF-8");
	static final byte[] BGEN_MAGIC = "bgen".getBytes(ENCODING);
	static final String FILE_SUFFIX = ".bgen";

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
    
    
	
	
	
	

	
	
	static abstract class RandomAccessStream extends FilterInputStream implements LocationAware
		{
		protected RandomAccessStream(InputStream in) {
			super(in);
			}
		public final long ftell()  throws IOException {
			return this.getPosition();	
			}
		public abstract boolean eof() throws IOException;
		public abstract void fseek(long n) throws IOException;
		}
	
	static class RandomAccessSeekableStream  extends RandomAccessStream {
		RandomAccessSeekableStream(SeekableStream in) {
			super(in instanceof SeekableBufferedStream?
					SeekableBufferedStream.class.cast(in):
					SeekableStreamFactory.getInstance().getBufferedStream(in)
					);
			}
		RandomAccessSeekableStream(String fname) throws IOException {
			this(SeekableStreamFactory.getInstance().getStreamFor(fname));
			}
		@Override
		public boolean eof()  throws IOException {
			return SeekableStream.class.cast(super.in).eof();
			}
		@Override
		public long getPosition() {
			try {
				return SeekableStream.class.cast(super.in).position();
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public void fseek(long position) throws IOException {
			SeekableStream.class.cast(super.in).seek(position);
			}
		}

	static class RandomAccessPositionalStream  extends RandomAccessStream {
		RandomAccessPositionalStream(InputStream in) {
			super(in instanceof PositionalBufferedStream?
					PositionalBufferedStream.class.cast(in):
					new PositionalBufferedStream(in)
					);
			}
		@Override
		public boolean eof() throws IOException {
			return PositionalBufferedStream.class.cast(super.in).isDone();
			}
		@Override
		public long getPosition() {
			return PositionalBufferedStream.class.cast(super.in).getPosition();
			}
		@Override
		public void fseek(long offset) throws IOException {
			final long position = getPosition();
			if(position>=offset) {
				skipNBytes(offset-position);
				}
			else
				{
				throw new IOException("random access seek operation not allowed when reading a stream");
				}
			}
		}
	
	
	public static class AlleleCombinations {
		private final Layout layout;
		private final List<String> alleles;
		private final int n_alleles;
		private final int m_ploidy;
		private final boolean phased;
		
		private static class FindByIndex implements Consumer<int[]> {
			final int wanted_index;
			int curr=0;
			int[] indexes = null;
			FindByIndex(final int wanted_index) {
				this.wanted_index = wanted_index;
				}
			@Override
			public void accept(int[] t) {
				if(curr==this.wanted_index) {
					this.indexes=Arrays.copyOf(t, t.length);
					}
				curr++;
				}
			}
		
		AlleleCombinations(Layout layout,int n_alleles,int m_ploidy,boolean phased) {
			this.layout=layout;
			this.alleles=null;
			this.n_alleles=  n_alleles;
			this.m_ploidy=m_ploidy;
			this.phased = phased;
			}
		public AlleleCombinations(Layout layout,final List<String> alleles,int m_ploidy,boolean phased) {
			this.layout=layout;
			this.alleles=new ArrayList<>(alleles);
			this.n_alleles = alleles.size();
			this.m_ploidy=m_ploidy;
			this.phased = phased;
			}
		 String getAllele(int idx) {
			 return this.alleles==null?
				 String.valueOf((char)('A'+idx)):
				 this.alleles.get(idx)
				 ;
		 	}
		
		 int calculateTotalCombinations() {
			 switch(layout) {
				 case LAYOUT_1: return 3;
				 case LAYOUT_2: return calculateTotalCombinationsForLayout2(phased,n_alleles,m_ploidy);
				 default: throw new IllegalArgumentException(); 
			 	}
			}
		
		 public int[] getAlleleIndexesByIndex(int i) {
			  if(layout.equals(Layout.LAYOUT_1)  || (layout.equals(Layout.LAYOUT_2) && !phased && n_alleles==2 && m_ploidy==2 )) {
				switch(i) {
					case 0: return new int[] {0,0};
					case 1: return new int[] {0,1};
					case 2: return new int[] {1,1};
					 default: throw new IllegalArgumentException(); 
					}
			  	}

			  
			  final int[] buffer=new int[this.m_ploidy];
			  FindByIndex andThen = new FindByIndex(i);
			  if(phased) {
		        	visitPhased(andThen, buffer, 0);
		        	}
		        else
		        	{
		        	visitUnphased(andThen, buffer, 0, 0);
		        	}
			  return andThen.indexes;
		  	}
		
		 
		  public List<int[]> getAllGenotypesIndexes() {
			  final List<int[]> container;
			if(layout.equals(Layout.LAYOUT_1)) {
				container = new ArrayList<>(3);
				container.add(new int[] {0,0});
				container.add(new int[] {0,1});
				container.add(new int[] {1,1});
			  } else {
			container = new ArrayList<>(calculateTotalCombinations());
        	final Consumer<int[]> anThen =A->container.add(Arrays.copyOf(A, A.length));
        	visit(anThen);
			  }
	        return container;
		  	}
		  
		public void debug() {
			final int[] tmp_idx=new int[]{0};
			final StringBuilder sb=new StringBuilder("Layout: "+layout+" Phased: "+phased+" ploidy: "+this.m_ploidy+" n-alleles:"+this.n_alleles+" {");
			if(layout.equals(Layout.LAYOUT_1)) {
				sb.append( " [0] : ").append(getAllele(0)).append(getAllele(0)).append(" | ");
				sb.append( " [1] : ").append(getAllele(0)).append(getAllele(1)).append(" | ");
				sb.append( " [1] : ").append(getAllele(1)).append(getAllele(1));
				}
			else
				{
				final Consumer<int[]> anThen =A->{
					sb.append( " ["+tmp_idx[0]+"] : ");
					for(int i=0;i< A.length;++i) {
						sb.append(getAllele(A[i]));
						}
					sb.append(" | ");
					tmp_idx[0]++;
					};
				visit(anThen);
				sb.append("}");
				}
	        LOG.debug(sb.toString());
			}
		
		private void visit(Consumer<int[]> consummer) {
			final int[] buffer=new int[this.m_ploidy];
	        if(phased) {
	        	visitPhased(consummer, buffer, 0);
	        	}
	        else
	        	{
	        	visitUnphased(consummer, buffer, 0, 0);
	        	}
			}
		  /*
		     alleles: 0,1,2
		     list:
		             0/0
		             0/1
		             0/2
		             1/1
		             1/2
		             2/2
		     */
	    private void visitUnphased(Consumer<int[]> anThen, int[] current, int start,int array_index) {	    	
	        if (array_index == m_ploidy) {
	        	anThen.accept(current);
	        	}
	        else
		        {
		        for (int i = start; i < n_alleles; i++) {
		            current[array_index]=i;
		            visitUnphased(anThen, current, i,array_index+1); // Allow repetition by not incrementing 'i'
		        	}
		        }
		    }
		  /*
	     alleles: 0,1,2
	     list:
	             0|0
	             0|1
	             0|2
	             1|0
	             1|1
	             1|2
	             2|0
	             2|1
	             2|2
	     */
	    private void visitPhased(Consumer<int[]> anThen, int[] current, int array_index) {
	        if (array_index == m_ploidy) {
	        	anThen.accept(current);
	        }
	        else {
	        for (int i = 0; i < n_alleles; i++) {
		            current[array_index] = i;
		            visitPhased(anThen, current, array_index + 1); // Move to the next depth
		        	}
		        }
	    	}
		}
	
}
