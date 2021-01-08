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

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.fastq.FastqPairedReaderFactory;
import com.github.lindenb.jvarkit.fastq.FastqPairedWriter;
import com.github.lindenb.jvarkit.fastq.FastqPairedWriterFactory;
import com.github.lindenb.jvarkit.fastq.FastqRecordPair;
import com.github.lindenb.jvarkit.fastq.FastqUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

/**

BEGIN_DOC

## Example

```bash


```

END_DOC

 */
@Program(name="fastqsplit",
	description="Split Fastq files into multiple files.",
	keywords="fastq",
	creationDate="20200327",
	modificationDate="20200327"
	)
public class FastqSplit extends Launcher
	{
	private static final Logger LOG = Logger.build(FastqSplit.class).make();
	private final static String TAG="__TOKEN__";

	@Parameter(names={"-o","--output"},description="Output file name. It MUST  end with a fastq suffix and MUST contain the word "+TAG,required=true)
	private String basename = null;
	@Parameter(names={"-s","-S","--splits"},description="number of splits. At most 'x' (pair-of)files will be generated. Or use option 'n' ")
	private int split_number = -1;
	@Parameter(names={"-n","-N","--count"},description="number of reads per file. will be generated. Or use option 's'.")
	private int per_file_number = -1;
	@Parameter(names={"-ii","--input-interleaved"},description="input is paired interleaved")
	private boolean input_is_interleaved = false;
	@Parameter(names={"-oi","--output-interleaved"},description="output is interleaved")
	private boolean output_is_interleaved = false;
	@Parameter(names={"-m","--manifest"},description="Optional manifest file")
	private Path manifestPath = null;
	@Parameter(names={"-md5","--md5"},description="write md5 file")
	private boolean write_md5=false;
	@Parameter(names={"-async","--async"},description="use async I/O")
	private boolean with_asynio=false;
	@Parameter(names={"-validate","--validate"},description="validate read names")
	private boolean validate_read_names=false;
	
	private FastqWriter openSingleWriter(int i,PrintWriter manifest) throws IOException {
		final FastqWriterFactory fqwf = new FastqWriterFactory();
		fqwf.setCreateMd5(this.write_md5);
		fqwf.setUseAsyncIo(this.with_asynio);
		
		final String tag = String.format("%09d.R0",(i+1));
		final String filename= this.basename.replace(TAG, tag);
		final File file = FastqUtils.validateFastqFilename(new File(filename));
		if(file.getParentFile()!=null) file.getParentFile().mkdirs();
		final FastqWriter w = fqwf.newWriter(file);
		manifest.println(file);
		return w;
		}
	
	private FastqPairedWriter openPairedWriter(int i,PrintWriter manifest) throws IOException {
	final FastqPairedWriterFactory fqpwf = new FastqPairedWriterFactory().
			setCreateMd5(write_md5).
			setAsyncIo(with_asynio);
	final FastqPairedWriter w;
	if(this.output_is_interleaved) {
		final String tag = String.format("%09d.R12",(i+1));
		final String filename= this.basename.replace(TAG, tag);
		final File file = new File(filename);
		if(file.getParentFile()!=null) file.getParentFile().mkdirs();
		w = fqpwf.open(file);
		manifest.println(file);
		}
	else
		{
		final File files[] = new File[2];
		for(int j=0;j< 2;++j) {
			final String tag = String.format("%09d.R%d",(i+1),(j+1));
			final String filename= this.basename.replace(TAG, tag);
			final File file = new File(filename);
			if(file.getParentFile()!=null) file.getParentFile().mkdirs();
			files[j]=file;
			}
		w = fqpwf.open(files[0],files[1]);
		manifest.print(files[0]);
		manifest.print("\t");
		manifest.println(files[1]);
		
		}

	return w;
	}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.per_file_number<1 && this.split_number < 1) {
			LOG.error("Option -n or -s are undefined");
			return -1;
			}
		if(this.per_file_number>0 && this.split_number > 0) {
			LOG.error("Both Options -n and -s are defined");
			return -1;
			}
		if(!this.basename.contains(TAG)) {
			LOG.error("basename doesn't contain "+TAG+": "+basename);
			return -1;
			}
		PrintWriter manifest=null;
		try
			{
			if(this.manifestPath==null) {
				manifest = new PrintWriter(new NullOuputStream());
			} else
			{
				manifest = super.openPathOrStdoutAsPrintWriter(this.manifestPath);
			}
			
			if(args.size()==2 || (this.input_is_interleaved && (args.isEmpty() || args.size()==1))) {
				final List<FastqPairedWriter> fastqWriters = new ArrayList<>();
				FastqPairedWriter previous=null;
				int count_files=0;
				long n = 0L;
				try(final CloseableIterator<FastqRecordPair> iter=   new FastqPairedReaderFactory().setValidateReadNames(this.validate_read_names).open(args)) {	
					while(iter.hasNext()) {
						final FastqRecordPair pair = iter.next();
						final FastqPairedWriter w ;
						/* split by number of reads per file */
						if(this.per_file_number>0) {
							if(previous==null || n%this.per_file_number==0) {
								if(previous!=null) previous.close();
								previous = openPairedWriter(count_files, manifest);
								count_files++;
								n=0L;
								}
							w = previous;
							}
						/* split by file */
						else 
							{
							final int idx = (int)(n%this.split_number);
							if(idx>= fastqWriters.size()) {
								w = openPairedWriter(idx, manifest);
								fastqWriters.add(w);
								}
							else
								{
								w = fastqWriters.get(idx);
								}
							}
						w.write(pair);
						n++;
						}
					}
				if(previous!=null) previous.close();
				for(final FastqPairedWriter w:  fastqWriters) w.close();
				}
			else if(args.isEmpty() || args.size()==1)
				{
				if(this.output_is_interleaved) {
					LOG.error("Cannot set output is interleaved if input is not paired.");
					return -1;
					}
				final List<FastqWriter> fastqWriters = new ArrayList<>();
				FastqWriter previous = null;
				long n=0L;
				int count_files = 0;
				final FastqReader iter;
				if(args.size()==1) {
					iter = new FastqReader(new File(args.get(0)));
					}
				else
					{
					iter = new FastqReader(IOUtils.openStreamForBufferedReader(stdin()));
					}	
				
				while(iter.hasNext()) {
					final FastqRecord rec = iter.next();
					final FastqWriter w ;
					
					/* split by number of reads per file */
					if(this.per_file_number>0) {
						if(previous==null || n%this.per_file_number==0) {
							if(previous!=null) previous.close();
							previous = this.openSingleWriter(count_files, manifest);
							count_files++;
							n=0L;
							}
						w = previous;
						}
					/* split by file */
					else 
						{
						final int idx = (int)(n%this.split_number);
						if( idx >= fastqWriters.size()) {
							w = this.openSingleWriter(idx, manifest);
							fastqWriters.add(w);
							}
						else
							{
							w = fastqWriters.get(idx);
							}
						}
					w.write(rec);
					n++;
					}
				iter.close();
				for(final FastqWriter w:  fastqWriters) w.close();
				}
			else
				{
				LOG.error("Illegal number of arguments.");
				return -1;
				}
			manifest.flush();
			manifest.close();
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(manifest);
			}
		}
	
	public static void main(final String[] args)
		{
		new FastqSplit().instanceMainWithExit(args);
		}

}
