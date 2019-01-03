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
package com.github.lindenb.jvarkit.io.tar;

import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigInteger;
import java.util.Arrays;

public class TarReader implements Closeable
	{
	private static final int BLOCK_SIZE=512;
	private static final byte zeroBlock[]=new byte[BLOCK_SIZE];
	{{{
		Arrays.fill(zeroBlock, (byte)0);
	}}}
	
	private static byte[] _malloc(int size)
		{
		return new byte[size];
		}
	/* Posix tar header*/
	private  class TARFileHeader
		{
		byte name[]=_malloc(100);
		byte mode[]=_malloc(8);
		byte uid[]=_malloc(8);
		byte gid[]=_malloc(8);
		byte size[]=_malloc(12);
		byte mtime[]=_malloc(12);
		byte checksum[]=_malloc(8);
		byte typeflag[]=_malloc(1);
		byte linkname[]=_malloc(100);
		byte magic[]=_malloc(6);
		byte version[]=_malloc(2);
		byte uname[]=_malloc(32);
		byte gname[]=_malloc(32);
		byte devmajor[]=_malloc(8);
		byte devminor[]=_malloc(8);
		byte prefix[]=_malloc(155);
		byte pad[]=_malloc(12);
		
		public String getName()
			{
			int length=0;
			while(length< name.length && name[length]!=0)
				{
				length++;
				}
			return new String(name,0,length);
			}
		
		public long getSize()
			{
			return parseOctalOrBinary(size,0,size.length);
			}
		}
	
	public  class Entry
		extends InputStream
		{
		private final long length;
		private long nRead=0L;
		private final TARFileHeader meta;
		private Entry(final TARFileHeader meta)
			{
			this.meta=meta;
			this.length=meta.getSize();
			}
		public long getOffset()
			{
			return nRead;
			}
		public long getLength()
			{
			return length;
			}
		@Override
		public int read() throws IOException {
			if(nRead>=length) return -1;
			int c=TarReader.this.in.read();
			if(c==-1) return -1;
			nRead++;
			if(nRead==length)
				{
				finish();
				}
			return c;
			}
		
		
		private void finish() throws IOException
			{
			while(nRead%BLOCK_SIZE!=0)
				{
				if(TarReader.this.in.read()==-1) throw new IOException();
				++nRead;
				}
			}
		
		public void skip() throws IOException
			{
			while(read()!=-1) ;
			}
		
		public String getName()
			{
			return meta.getName();
			}
		
		@Override
		public void close()
			{
			}
		@Override
		public String toString() {
			return getName()+" "+nRead+"/"+length;
			}
		}
	
	
	private InputStream in=null;
	public TarReader(InputStream in) throws IOException
		{
		this.in=in;
		}
	static int copy(byte array[],byte field[],int off)
		{
		System.arraycopy(array, off,field,0,field.length);
		return off+field.length;
		}
	private TARFileHeader readTARFileHeader() throws IOException
		{
		byte array[]=new byte[BLOCK_SIZE];
		int nRead=0;
		while(nRead< BLOCK_SIZE)
			{
			int x=in.read(array,nRead,BLOCK_SIZE-nRead);
				
			if(x<0)
				{
				throw new IOException("cannot read bytes "+nRead+" / "+array.length);
				}
			nRead+=x;
			}
		if(Arrays.equals(array, zeroBlock)) return null;
		TARFileHeader h=new TARFileHeader();
		int off=0;
		off=copy(array,h.name,off);
		off=copy(array,h.mode,off);
		off=copy(array,h.uid,off);
		off=copy(array,h.gid,off);
		off=copy(array,h.size,off);
		off=copy(array,h.mtime,off);
		off=copy(array,h.checksum,off);
		off=copy(array,h.typeflag,off);
		off=copy(array,h.linkname,off);
		off=copy(array,h.magic,off);
		off=copy(array,h.version,off);
		off=copy(array,h.uname,off);
		off=copy(array,h.gname,off);
		off=copy(array,h.devmajor,off);
		off=copy(array,h.devminor,off);
		off=copy(array,h.prefix,off);
		off=copy(array,h.pad,off);
		if(off!=BLOCK_SIZE) throw new IOException();
		return h;
		}
	@Override
	public void close() throws IOException {
		if(in!=null) in.close();
		in=null;
		}
	
	public Entry next() throws IOException
		{
		TARFileHeader h= readTARFileHeader();
		if(h==null) return null;
		return new Entry(h);
		}	
/* 
-rw-r--r-- tm/tm      11255269 2014-01-23 14:20 citations.dmp
-rw-r--r-- tm/tm       2768900 2014-01-23 14:20 delnodes.dmp
-rw-r--r-- tm/tm           419 2014-01-23 14:20 division.dmp
-rw-r--r-- tm/tm         11559 2014-01-23 14:20 gc.prt
-rw-r--r-- tm/tm          3566 2014-01-23 14:20 gencode.dmp
-rw-r--r-- tm/tm        554847 2014-01-23 14:20 merged.dmp
-rw-r--r-- tm/tm      96596932 2014-01-23 14:20 names.dmp
 */
	
    private static long parseOctalOrBinary(final byte[] buffer, final int offset,
            final int length) {
			
			if ((buffer[offset] & 0x80) == 0) {
			return parseOctal(buffer, offset, length);
			}
			final boolean negative = buffer[offset] == (byte) 0xff;
			if (length < 9) {
			return parseBinaryLong(buffer, offset, length, negative);
			}
			return parseBinaryBigInteger(buffer, offset, length, negative);
			}
    
    private static long parseOctal(final byte[] buffer, final int offset, final int length) {
        long    result = 0;
        int     end = offset + length;
        int     start = offset;

        if (length < 2){
            throw new IllegalArgumentException("Length "+length+" must be at least 2");
        }

        if (buffer[start] == 0) {
            return 0L;
        }

        // Skip leading spaces
        while (start < end){
            if (buffer[start] == ' '){
                start++;
            } else {
                break;
            }
        }

        // Must have trailing NUL or space
        byte trailer;
        trailer = buffer[end-1];
        if (trailer == 0 || trailer == ' '){
            end--;
        } else {
            throw new IllegalArgumentException("");
        }
        // May have additional NUL or space
        trailer = buffer[end-1];
        if (trailer == 0 || trailer == ' '){
            end--;
        }

        for ( ;start < end; start++) {
            final byte currentByte = buffer[start];
            // CheckStyle:MagicNumber OFF
            if (currentByte < '0' || currentByte > '7'){
                throw new IllegalArgumentException("");
            }
            result = (result << 3) + (currentByte - '0'); // convert from ASCII
            // CheckStyle:MagicNumber ON
        }

        return result;
    }

    
    private static long parseBinaryLong(final byte[] buffer, final int offset,
                final int length,
                final boolean negative)
			  {
			if (length >= 9) {
			throw new IllegalArgumentException("At offset " + offset + ", "
			                       + length + " byte binary number"
			                       + " exceeds maximum signed long"
			                       + " value");
			}
			long val = 0;
			for (int i = 1; i < length; i++) {
			val = (val << 8) + (buffer[offset + i] & 0xff);
			}
			if (negative) {
			// 2's complement
			val--;
			val ^= ((long) Math.pow(2, (length - 1) * 8) - 1);
			}
			return negative ? -val : val;
			  }
        private static long parseBinaryBigInteger(final byte[] buffer,
                final int offset,
                final int length,
                final boolean negative) {
				byte[] remainder = new byte[length - 1];
				System.arraycopy(buffer, offset + 1, remainder, 0, length - 1);
				BigInteger val = new BigInteger(remainder);
				if (negative) {
				// 2's complement
				val = val.add(BigInteger.valueOf(-1)).not();
				}
				if (val.bitLength() > 63) {
				throw new IllegalArgumentException("At offset " + offset + ", "
				                 + length + " byte binary number"
				                 + " exceeds maximum signed long"
				                 + " value");
				}
				return negative ? -val.longValue() : val.longValue();
				}

	}
