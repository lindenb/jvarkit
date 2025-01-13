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
package com.github.lindenb.jvarkit.fastq;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;

import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

import htsjdk.samtools.fastq.FastqRecord;

public class FastqRecordCodec  extends AbstractDataCodec<FastqRecord> {
	
	private static String notNull(final String s)
		{
		return s==null?"":s;
		}	

	
	@Override
	public FastqRecord decode(final DataInputStream dis) throws IOException {
		final String seqHeader;
		try {
			seqHeader = dis.readUTF();
			}
		catch(EOFException err) {
			return null;
			}
		final String qualHeader=dis.readUTF();
		final int len = dis.readInt();
		final byte array[] = new byte[len*2];
		dis.readFully(array);
		return new FastqRecord(seqHeader,
				new String(array,0,len),
				qualHeader ,
				new String(array,len,len)
				);
		}

	@Override
	public void encode(final DataOutputStream dos, final FastqRecord r) throws IOException {
		dos.writeUTF(notNull(r.getReadName()));
		dos.writeUTF(notNull(r.getBaseQualityHeader()));
		dos.writeInt(r.getReadLength());
		dos.write(r.getReadBases());
		dos.write(r.getBaseQualities());
		}

	@Override
	public AbstractDataCodec<FastqRecord> clone() {
		return new FastqRecordCodec();
	}
	
}
