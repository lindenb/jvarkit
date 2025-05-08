package com.github.lindenb.jvarkit.bgen;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.charset.Charset;
import java.util.function.DoubleSupplier;

import com.github.lindenb.jvarkit.io.IOUtils;

import htsjdk.samtools.Defaults;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.tribble.TribbleException;

public class BGenUtils {
	BGenUtils() {
		
	}
	
	static final Charset ENCODING=Charset.forName("UTF-8");
	static final byte[] BGEN_MAGIC = "bgen".getBytes(ENCODING);
	static final String FILE_SUFFIX = ".bgen";
	static final int USHRT_MAX= 65_635;

    /** convert a long to int, throws a tribble exception if the number is greater than Integer.MAX_VALUE */
    static int longToUnsignedInt(final long n) {
        if(n>Integer.MAX_VALUE || n<0L) throw new TribbleException("Cannot convert "+n+"L to int");
        return (int)n;
    	}
    /** read the length of a string as unsigned short and read the string */
    static String readStringUInt16(final BinaryCodec binaryCodec) throws IOException {
    	return readString(binaryCodec, binaryCodec.readUShort());
    	}
    /** read the length of a string as unsigned short and read the string */
    public static String readStringUInt32(final BinaryCodec binaryCodec) throws IOException {
    	return readString(binaryCodec, longToUnsignedInt(binaryCodec.readUInt()));
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

   	/** ByteArrayOutputStream but we expose the internal buffer */
    protected static class ByteBuffer extends ByteArrayOutputStream {
		byte[] getBuffer() { return super.buf;}
		int getCount() { return super.count;}
		public void copyTo(InputStream in,int nbytes) throws IOException {
			IOUtils.copyTo(in,this,nbytes);
			}
		ByteArrayInputStream toByteArrayInputStream() {
			return new ByteArrayInputStream(getBuffer(), 0, getCount());
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
				t=t*2 + c;
				--n;
				}
			return t/max;
			}
		}
}
