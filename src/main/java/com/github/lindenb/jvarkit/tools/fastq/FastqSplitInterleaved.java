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
package com.github.lindenb.jvarkit.tools.fastq;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.fastq.FastqPairedReaderFactory;
import com.github.lindenb.jvarkit.fastq.FastqPairedWriter;
import com.github.lindenb.jvarkit.fastq.FastqPairedWriterFactory;
import com.github.lindenb.jvarkit.fastq.FastqRecordPair;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

/**

BEGIN_DOC

## Example

```bash
$ curl -sk "https://raw.githubusercontent.com/bigdatagenomics/adam/fff8ae259e8f6958eefd8de9a3ec39d33392fb21/adam-core/src/test/resources/interleaved_fastq_sample1.fq" |\
	java -jar dist/fastqsplitinterleaved.jar  -R1 R1.fastq.gz -R2 R1.fastq.gz 
```

END_DOC

 */
@Program(name="fastqsplitinterleaved",
	description="Split interleaved Fastq files.",
	keywords="fastq",
	modificationDate="20200327"
	)
public class FastqSplitInterleaved extends Launcher
	{
	private static final Logger LOG = Logger.build(FastqSplitInterleaved.class).make();
	@Parameter(names={"-a","-R1","-F"},description="R1 file",required=true)
	private Path fileA=null;
	@Parameter(names={"-b","-R2","-R"},description="R2 file",required=true)
	private Path fileB=null;
	@Parameter(names={"-md5","--md5"},description="write md5 file")
	private boolean write_md5=false;
	@Parameter(names={"-async","--async"},description="use async I/O")
	private boolean with_asynio=false;
	@Parameter(names={"-validate","--validate"},description="validate read names")
	private boolean validate_read_names=false;

	
	private FastqSplitInterleaved() {
		}
	

	@Override
	public int doWork(final List<String> args) {
		if(this.fileA.equals(this.fileB)) {
			LOG.error("R1 file==R2.file.");
			return -1;
			}
		CloseableIterator<FastqRecordPair> iter1=null;
		FastqPairedWriter pairedWriter = null;
		try
			{
			final String input = oneFileOrNull(args);
			final FastqPairedReaderFactory fqprf = new FastqPairedReaderFactory().setValidateReadNames(this.validate_read_names);
			if(input==null)
				{
				iter1 = fqprf.open(stdin());	
				}
			else
				{
				iter1 = fqprf.open(Paths.get(input));				
				}
			
			final FastqPairedWriterFactory fqwf = new FastqPairedWriterFactory().
					setCreateMd5(this.write_md5).
					setAsyncIo(this.with_asynio);
			pairedWriter = fqwf.open(fileA, fileB);
			
			while(iter1.hasNext())
				{
				pairedWriter.write(iter1.next());				
				}
			iter1.close();
			pairedWriter.close();
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter1);
			CloserUtil.close(pairedWriter);
			}
		}
	
	public static void main(final String[] args)
		{
		new FastqSplitInterleaved().instanceMainWithExit(args);
		}
	}
