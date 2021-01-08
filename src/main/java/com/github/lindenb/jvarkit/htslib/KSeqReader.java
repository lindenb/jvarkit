/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.htslib;

import java.io.IOException;

import com.github.lindenb.jvarkit.jni.CPtr;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloseableIterator;

public class KSeqReader extends CPtr implements CloseableIterator<FastqRecord>{
	private final String filename;
	private FastqRecord curr = null;
	private final String array4[] = new String[4];
	public KSeqReader(long ptr, final String filename) {
		super(ptr);
		this.filename = filename;
	}

	@Override
	public boolean hasNext() {
		if(isNull()) return false;
		if(this.curr!=null) return true;
		
		final int rez = HtsLib.kseq_read4(getPtr(),this.array4);
		if(rez<0) return false;
		this.curr = new FastqRecord(
				this.array4[0],
				this.array4[1],
				this.array4[2],
				this.array4[3]
				); 
		return true;
		}
	@Override
	public FastqRecord next() {
		if(!hasNext()) throw new IllegalStateException("no more read in "+this.filename);
		FastqRecord rec = this.curr;
		this.curr=null;
		return rec;
		}
	
	@Override
	public void close() {
		if(!isNull()) HtsLib.kseq_destroy(getPtr());
		setNull();
		this.curr=null;
		}
	
	public static KSeqReader open(final String filename) throws IOException {
		final long ptr = HtsLib.kseq_init_file(filename);
		if(ptr==0L) throw new IOException("Cannot open kseq file \""+filename+"\"");
		return new KSeqReader(ptr,filename);
		}
}
