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

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.fastq.FastqPairedReaderFactory;
import com.github.lindenb.jvarkit.fastq.FastqRecordPair;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.util.CloseableIterator;

public abstract class FastqLauncher extends Launcher {
private static final Logger LOG = Logger.of(FastqLauncher.class);
@Parameter(names={"--paired"},description="assume input is paired end: we expect two files, or the input is assumed interleaved fastq.")
boolean paired_end = false;


protected Logger getLogger() {
	return LOG;
	}

protected int beforeFastq() {
	return 0;
	}
protected void afterFastq() {
	}

protected int runPairedEnd(final CloseableIterator<FastqRecordPair> iter) throws IOException {
	throw new UnsupportedOperationException("Paired-end not implemented");
	}

protected int runSingleEnd(final FastqReader fqr) throws IOException {
	throw new UnsupportedOperationException("Single-end not implemented");
}

/** does this software supported paired-end n_ends==2 or singled_end n_ends==1 */
protected boolean isSupportingXXXEnd(int n_ends) {
	return true;
	}

@Override
public int doWork(final List<String> args)
	{
	try {
	if(beforeFastq()!=0) {
		getLogger().error("initialization failed.");
		return 1;
	}
	final int ret;
	if(paired_end || args.size()==2) {
		if(!paired_end && args.size()==2) {
			getLogger().warn("two files for input. Assuming paired-end input");
			}
		if(!isSupportingXXXEnd(2)) {
			getLogger().error("This software doesn't support paired-end input");
			return -1;
			}
		final FastqPairedReaderFactory fqpr = new FastqPairedReaderFactory();
		try(CloseableIterator<FastqRecordPair> iter=fqpr.open(args)) {
			ret = runPairedEnd(iter);
			}
		}
	else
		{
		if(!isSupportingXXXEnd(1)) {
			getLogger().error("This software doesn't support single-end input");
			return -1;
			}
		final String input = oneFileOrNull(args);
		try(FastqReader fqr= (input==null?
				new FastqReader(IOUtils.openStdinForBufferedReader()):
				new FastqReader(new File(input))
				)){
				ret = runSingleEnd(fqr);
				}
			}
		return ret;
		} catch(final Throwable err) {
		
		
		getLogger().error(err);
		return -1;
		}
	finally {
		try { afterFastq();} catch(final Throwable err2) {
			getLogger().warn(err2);
			}
		}
	}

}
