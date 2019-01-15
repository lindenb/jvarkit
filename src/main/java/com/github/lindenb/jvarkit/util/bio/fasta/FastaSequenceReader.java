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


History:
* 2016 creation

*/
package com.github.lindenb.jvarkit.util.bio.fasta;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PushbackReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;

public class FastaSequenceReader {
private int sequenceCapacity=100;
private BiFunction<String,byte[],FastaSequence> fastaSequenceCreator = new BiFunction<String, byte[], FastaSequence>() {
	@Override
	public FastaSequence apply(final String name,final  byte[] seq) {
		return new SequenceImpl(name, seq);
	}
};


public FastaSequenceReader() {
	}

public void setSequenceCapacity(final int sequenceCapacity) {
	this.sequenceCapacity = sequenceCapacity;
}

public int getSequenceCapacity() {
	return sequenceCapacity;
	}



public void setFastaSequenceCreator(BiFunction<String, byte[], FastaSequence> fastaSequenceCreator) {
	this.fastaSequenceCreator = fastaSequenceCreator;
	}	

public CloseableIterator<FastaSequence> iterator(final Reader r) throws IOException {
	return new MyIterator(r);
	}

public CloseableIterator<FastaSequence> iterator(final File file) throws IOException {
	return iterator(IOUtils.openFileForReader(file));
	}

public Iterable<FastaSequence> getSequencesIn(final File file)  {
	return new FastaIterable(file);
	}


/** read one sequence in the fasta file. There must be one AND ONLY one sequence */
public FastaSequence readOne(final File file) throws IOException {
	CloseableIterator<FastaSequence> r= null;
	try
		{
		r =	this.iterator(file);
		if(!r.hasNext()) {
			r.close();r=null;
			throw new IOException("Expected one sequence in "+file+" but got none");
			}
		final FastaSequence seq = r.next();
		if(r.hasNext()) {
			r.close();r=null;
			throw new IOException("Expected only one sequence in "+file+" but got none after "+seq.getName() );
			}
		r.close();r=null;
		return seq;
		}
	finally
		{
		CloserUtil.close(r);
		}
	}

public List<FastaSequence> readAll(final File file) throws IOException {
	final List<FastaSequence> seqs = new ArrayList<>(); 
	final CloseableIterator<FastaSequence> r= iterator(file);
	while(r.hasNext())  seqs.add(r.next());
	r.close();
	return seqs;
	}

protected FastaSequence createFastaSequence(final String name,byte seq[]) {
	return new SequenceImpl(name, seq);
	}

protected FastaSequence read(final PushbackReader reader) throws IOException {
	boolean at_begin=true;
	StringBuilder name=null;
	ByteArrayOutputStream sequence=null;
	try {
		int c;
		while((c=reader.read())!=-1) {
			if(at_begin && c=='>') {
				if(name!=null) {
					reader.unread(c);
					return this.fastaSequenceCreator.apply(
							name.toString(),
							sequence.toByteArray()
							);
					}
				name = new StringBuilder();
				sequence =new ByteArrayOutputStream(this.sequenceCapacity);
				/* consume header */
				while((c=reader.read())!=-1 && c!='\n') {
				name.append((char)c);	
				}
				at_begin = true;
			} else if(Character.isWhitespace(c)) {
				at_begin = (c=='\n');
			} else if(sequence==null) {
				throw new IOException("Illegal character "+(char)c);
			} else
			{
				sequence.write(c);	
				while((c=reader.read())!=-1 && c!='\n') {
					sequence.write(c);	
				}
				at_begin=true;
			}
		}
		/* eof met */
		if(name!=null) {
			return this.fastaSequenceCreator.apply(
					name.toString(),
					sequence.toByteArray()
					);
			}
		return null;
	} catch (final IOException e) {
		throw new RuntimeIOException(e);
		}}

private class MyIterator extends AbstractIterator<FastaSequence>
	implements CloseableIterator<FastaSequence>
	{
	PushbackReader r;
	MyIterator(final Reader r) {
		if(r==null) throw new RuntimeIOException("reader is null");
		this.r=new PushbackReader(r);
		}

	
	@Override
	protected FastaSequence advance() {
		if(r==null) return null;
		try {
			final FastaSequence s= FastaSequenceReader.this.read(this.r);
			if(s==null) close();
			return s;
		} catch (IOException e) {
			throw new RuntimeIOException();
		}
		}
	
	@Override
	public void close() {
		CloserUtil.close(r);
		r=null;
		}
	}

private static class SequenceImpl
	extends AbstractCharSequence
	implements FastaSequence {
	final String name;
	final byte seq[];
	SequenceImpl(final String name,final byte seq[]) {
		this.name=name;
		this.seq=seq;
	}
	@Override
	public String getName() {
		return name;
		}
	@Override
	public int length() {
		return seq.length;
		}
	@Override
	public char charAt(final int index) {
		return (char)seq[index];
		}
	}

	public class FastaIterable implements Iterable<FastaSequence> {
		final File fastaFile;
		FastaIterable(final File fastaFile) {
			this.fastaFile = fastaFile;
		}
		@Override
		public CloseableIterator<FastaSequence> iterator() {
			try {
				return FastaSequenceReader.this.iterator(this.fastaFile);
			} catch (final IOException e) {
				throw new RuntimeIOException(e);
			}
		}
	}

}
