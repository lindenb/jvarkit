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
package com.github.lindenb.jvarkit.jcommander;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.fastq.FastqPairedWriter;
import com.github.lindenb.jvarkit.fastq.FastqPairedWriterFactory;
import com.github.lindenb.jvarkit.fastq.FastqRecordPair;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.CloseableIterator;

public abstract class OnePassFastqLauncher extends FastqLauncher {
private static final Logger LOG = Logger.of(OnePassFastqLauncher.class);
@Parameter(names={"-o","--out","-R1"},description="Output file for R1 fastq record or interleaved output."+OPT_OUPUT_FILE_OR_STDOUT)
private File outputFile1 = null;
@Parameter(names={"-R2"},description="Output file for R2 fastq record. If input is paired but R2 is omitted, output will be interleaved.")
private File outputFile2 = null;
@Parameter(names={"-md5","--md5"},description="write md5 file")
private boolean write_md5=false;

@Override
protected Logger getLogger() {
	return LOG;
	}

protected Predicate<FastqRecordPair> createPredicateForFastqRecordPair() {
	throw new UnsupportedOperationException();
	}
protected int runPairedEnd(CloseableIterator<FastqRecordPair> iter,FastqPairedWriter fws) throws IOException {
	final  Predicate<FastqRecordPair> predicate = this.createPredicateForFastqRecordPair();
	while(iter.hasNext()) {
		final FastqRecordPair pair= iter.next();
		if(predicate.test(pair)) {
			fws.write(pair);
			}
		}
	return 0;
	}

protected Predicate<FastqRecord> createPredicateForFastqRecord() {
	throw new UnsupportedOperationException();
	}
protected int runSingleEnd(FastqReader iter,FastqWriter fws) throws IOException {
	final  Predicate<FastqRecord> predicate = this.createPredicateForFastqRecord();
	while(iter.hasNext()) {
		final FastqRecord rec= iter.next();
		if(predicate.test(rec)) {
			fws.write(rec);
			}
		}
	return 0;
	}

@Override
protected int runPairedEnd(final CloseableIterator<FastqRecordPair> iter) throws IOException {
	int ret = 0;
	FastqPairedWriter fws = null;
	try {
		final FastqPairedWriterFactory fpwf = new FastqPairedWriterFactory();
		fpwf.setCreateMd5(this.write_md5);
		
		if(outputFile1!=null && outputFile2!=null) {
			fws = fpwf.open(outputFile1,outputFile2);
			}
		else if(outputFile1!=null && outputFile2==null) {
			fws = fpwf.open(outputFile1);
			}
		else if(outputFile1==null && outputFile2==null) 
			{
			fws = fpwf.open(new PrintStream(new BufferedOutputStream(stdout())));
			}
		else
			{
			getLogger().error("bad output declaration.");
			return -1;
			}
		ret = runPairedEnd(iter, fws);
		fws.close();
		return ret;
		}
	catch(final Throwable err) {
		getLogger().error(err);
		return -1;
		}
	finally {
		if(fws!=null) fws.close();
		}
	}

@Override
protected int runSingleEnd(final FastqReader fqr) throws IOException {
if(this.outputFile1==null) {
	try(FastqWriter fw = new BasicFastqWriter(new PrintStream(new BufferedOutputStream(stdout())))){
		return runSingleEnd(fqr, fw);
		}
	} else {
		final FastqWriterFactory fqwf = new FastqWriterFactory();
		fqwf.setCreateMd5(this.write_md5);
		try(FastqWriter fw = fqwf.newWriter(this.outputFile1)) {
			return runSingleEnd(fqr, fw);
			}
		}
}

}
