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

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Path;

import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;

public class FastqPairedWriterFactory {
private boolean createMd5 = false;
private boolean useAsyncIo = false;

public FastqPairedWriterFactory setCreateMd5(boolean createMd5) {
	this.createMd5 = createMd5;
	return this;
	}
public boolean isCreateMd5() {
	return createMd5;
	}

public FastqPairedWriterFactory setAsyncIo(boolean useAsyncIo) {
	this.useAsyncIo = useAsyncIo;
	return this;
	}
public boolean isAsyncIo() {
	return useAsyncIo;
	}

/** write interleaved reads into two paths */
public FastqPairedWriter open(final Path f1,final Path f2) throws IOException {
	return open(f1.toFile(),f2.toFile());
	}

private FastqWriterFactory newFastqWriterFactory() {
	final FastqWriterFactory qfw = new FastqWriterFactory();
	qfw.setCreateMd5(isCreateMd5());
	qfw.setUseAsyncIo(isAsyncIo());
	return qfw;
}

/** write interleaved reads into two files */
public FastqPairedWriter open(final File f1,final File f2) throws IOException {
	final FastqWriterFactory qfw = newFastqWriterFactory();
	final FastqWriter w1 = qfw.newWriter(FastqUtils.validateFastqFilename(f1));
	final FastqWriter w2 = qfw.newWriter(FastqUtils.validateFastqFilename(f2));
	return new TwoWriter(w1,w2);
	}


/** write interleaved reads in one path */
public FastqPairedWriter open(final Path path) throws IOException {
	return open(path.toFile());
	}

/** write interleaved reads in one file */
public FastqPairedWriter open(final File fqFile) throws IOException {
	final FastqWriterFactory qfw = newFastqWriterFactory();
	final FastqWriter w = qfw.newWriter(fqFile);
	return new OneWriter(w);
	}

/** write interleaved reads to stream  */
public FastqPairedWriter open(final PrintStream w) throws IOException {
	return new OneWriter(new BasicFastqWriter(w));
	}



private static abstract class AbstractPairWriter implements FastqPairedWriter {
	final protected FastqWriter fq1w;
	final protected FastqWriter fq2w;
	protected AbstractPairWriter(final FastqWriter fqReader1,final FastqWriter fqReader2) {
		this.fq1w = fqReader1;
		this.fq2w = fqReader2;
		}
	@Override
	public void write(final FastqRecord r1, final FastqRecord r2) {
		this.fq1w.write(r1);
		this.fq2w.write(r2);
		}
	
	}

static class OneWriter extends AbstractPairWriter {
	OneWriter(final FastqWriter fqw) {
		super(fqw,fqw);
		}
	
	@Override
	public void close() {
		super.fq1w.close();
		}
	}

private static class TwoWriter extends AbstractPairWriter {
	TwoWriter(final FastqWriter w1,final FastqWriter w2) {
		super(w1,w2);
		}
	@Override
	public void close() {
		super.fq1w.close();
		super.fq2w.close();
		}
	}

}
