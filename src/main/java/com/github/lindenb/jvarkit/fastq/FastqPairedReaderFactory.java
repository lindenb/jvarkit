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
package com.github.lindenb.jvarkit.fastq;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;

public class FastqPairedReaderFactory {
	private boolean validate_read_names = false;
	
	public FastqPairedReaderFactory setValidateReadNames(boolean validate_read_names) {
		this.validate_read_names = validate_read_names;
		return this;
		}
	
	public boolean isValidateReadNames() {
		return validate_read_names;
		}
	
	/** 
	if args is empty, read from interleaved stdin
	if args.size==1 read from one interleaved file
	if args.size==2 read from two files */
	public CloseableIterator<FastqRecordPair> open(final List<String> args) throws IOException {
	switch(args.size()) {
		case 0: return open(System.in);
		case 1: return open(Paths.get(args.get(0)));
		case 2: return open(Paths.get(args.get(0)),Paths.get(args.get(1)));
		default: throw new IllegalArgumentException("exected 0,1 or 2 arguments but got "+args);
		}
	}
	
	/** read reads in two files */
	public CloseableIterator<FastqRecordPair> open(final Path fq1,final Path fq2) throws IOException {
		return open(fq1.toFile(),fq2.toFile());
	}
	
	/** read reads in two files */
	public CloseableIterator<FastqRecordPair> open(final File fq1,final File fq2) throws IOException {
		if(fq1.equals(fq2)) throw new IllegalArgumentException("same file "+fq1+"=="+fq2);
		return new InterleavedTwoReader(
				 new FastqReader(fq1),
				 new FastqReader(fq2)
				);
		}

/** read interleaved reads in one path */
public CloseableIterator<FastqRecordPair> open(final Path path) throws IOException {
	return open(path.toFile());
	}

	
/** read interleaved reads in one file */
@SuppressWarnings("resource")
public CloseableIterator<FastqRecordPair> open(final File fqFile) throws IOException {
	return new InterleavedOneReader( new FastqReader(fqFile)).setValidateReadName(this.isValidateReadNames());
	}

/** read interleaved reads from input stream */
public CloseableIterator<FastqRecordPair> open(final InputStream is) throws IOException {
	return open(new BufferedReader(new InputStreamReader(is)));
	}
/** read interleaved reads from buffered reader  */
@SuppressWarnings("resource")
public CloseableIterator<FastqRecordPair> open(final BufferedReader br) throws IOException {
	return new InterleavedOneReader( new FastqReader(br)).setValidateReadName(this.isValidateReadNames());
	}



private static abstract class AbstractInterleavedPairReader extends AbstractIterator<FastqRecordPair> 
implements CloseableIterator<FastqRecordPair> {
	final protected FastqReader fq1Reader;
	final protected FastqReader fq2Reader;
	protected boolean validate_read_name =false;
	protected AbstractInterleavedPairReader(final FastqReader fqReader1,final FastqReader fqReader2) {
		this.fq1Reader = fqReader1;
		this.fq2Reader = fqReader2;
		}
	
	AbstractInterleavedPairReader setValidateReadName(boolean b) {
		this.validate_read_name = b;
		return this;
	}
	
@Override
protected FastqRecordPair advance() {
	if(!this.fq1Reader.hasNext())  {
		if(this.fq2Reader.hasNext()) {
			final FastqRecord r2 = this.fq2Reader.next();
			close();
			throw new RuntimeIOException(
					"Expected Read in file "+this.fq1Reader.getFile()
					+ " after reading read "+r2+" in file "+this.fq2Reader.getFile());
			
			}
		return null;
		}
	final FastqRecord r1 = this.fq1Reader.next();
	if(!this.fq2Reader.hasNext()) {
		close();
		throw new RuntimeIOException(
				"Expected Read in file "+this.fq2Reader.getFile()
				+ " after reading read "+r1+" in file "+this.fq1Reader.getFile());
		}
	final FastqRecord r2 = this.fq2Reader.next();
	
	final FastqRecordPair pair =  FastqRecordPair.of(r1,r2);
	
	if(this.validate_read_name) pair.validateName();
	
	return pair;
	}
}


private static class InterleavedOneReader extends AbstractInterleavedPairReader {
	InterleavedOneReader(final FastqReader fqReader) {
		super(fqReader,fqReader);
		}
	@Override
	public void close() {
		super.fq1Reader.close();
		}
	}

private static class InterleavedTwoReader extends AbstractInterleavedPairReader {
	InterleavedTwoReader(final FastqReader r1,final FastqReader r2) {
		super(r1,r2);
		}
	@Override
	public void close() {
		super.fq1Reader.close();
		super.fq2Reader.close();
		}
	}

}
