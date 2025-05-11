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


History:
* 2015 : moving to knime
* 2014 : creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.List;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.tribble.readers.LineIterator;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
/**

BEGIN_DOC

## Deprecated

just use htslib/tabix...

## Example

```bash
  java -jar dist/bedindextabix.jar -s -o out.bed.gz input.bed

```

END_DOC

 */
@Program(
	name="bedindextabix",
	description="Index and sort a Bed on the fly with Tabix (deprecated).",
	keywords={"bed","tabix"},
	creationDate = "20150708",
	modificationDate = "20240724",
	jvarkit_amalgamion = true
	)
public class BedIndexTabix
	extends Launcher
	{
	private static final Logger LOG = Logger.of(BedIndexTabix.class);


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT,required=true)
	private Path outputFile = null;

	@Parameter(names={"-s","--sort"},description="sort BED prior to saving")
	private boolean sort=false;
	
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	private final BedLineCodec bedCodec=new BedLineCodec();
	
	
	
	private static class BedDataCodec extends  AbstractDataCodec<String>
			{
			@Override
			public String decode(DataInputStream dis)
					throws IOException
				{
				try {
					String s= AbstractDataCodec.readString(dis);
					return s;
				} catch (Exception e) {
					return null;
					}
				}
			@Override
			public void encode(DataOutputStream dos, String s) throws IOException {
				AbstractDataCodec.writeString(dos, s);
					}
			@Override
			public BedDataCodec clone() {
					return new BedDataCodec();
					}
			};

	
	public BedIndexTabix()
		{
		}
	
	
	
	protected void run(LineIterator in) throws IOException
		{
		int bedLineCount=0;
		
		final Path tbi = outputFile.getParent().resolve(outputFile.getFileName().toString()+FileExtensions.TABIX_INDEX);
		
		final Comparator<String> comparator=new Comparator<String>	()
			{
			@Override
			public int compare(String o1, String o2)
				{
				BedLine bed1 = bedCodec.decode(o1);
				BedLine bed2 = bedCodec.decode(o2);
				int i= bed1.getContig().compareTo(bed2.getContig());
				if(i!=0) return i;
				i= bed1.getStart() - bed2.getStart();
				if(i!=0) return i;
				i= bed1.getEnd() - bed2.getEnd();
				if(i!=0) return i;
				return o1.compareTo(o2);
				}
			};

		try {
			final TabixIndexCreator indexCreator=new TabixIndexCreator(TabixFormat.BED);
			LOG.info("Opening"+outputFile);
			try(BlockCompressedOutputStream writer = new BlockCompressedOutputStream(this.outputFile,BlockCompressedOutputStream.getDefaultCompressionLevel(),BlockCompressedOutputStream.getDefaultDeflaterFactory())) {
				final StringBuilder header=new StringBuilder();
				while(in.hasNext())
					{
					String h=in.peek();
					if(!BedLine.isBedHeader(h)) break;
					header.append(in.next()).append('\n');
					}
				//write header
				if(header.length()>0)
					{
					LOG.info("Writing header");
					writer.write(header.toString().getBytes());
					}
				
				if(this.sort)
					{
					LOG.info("Sorting");
					
					SortingCollection<String>  sorter =  SortingCollection.newInstance(
			                        String.class,
			                        new BedDataCodec(),
			                        comparator,
			                        this.writingSortingCollection.getMaxRecordsInRam(),
			                        this.writingSortingCollection.getTmpPaths()
			                        );
					while(in.hasNext())
						{
						String line = in.next();
						BedLine bed = bedCodec.decode(line);
						if(bed==null) continue;
						sorter.add(line);
						}
					sorter.doneAdding();
					sorter.setDestructiveIteration(true);
					try(CloseableIterator<String> iter= sorter.iterator()) {
						long filePosition= writer.getFilePointer();
						while(iter.hasNext())
							{
							String line = iter.next();
							BedLine bed = this.bedCodec.decode(line);
							writer.write(line.getBytes());
							writer.write('\n');
							indexCreator.addFeature(bed, filePosition);
							filePosition = writer.getFilePointer();
							}
						}
					sorter.cleanup();
					}
				else
					{
					long filePosition= writer.getFilePointer();
					while(in.hasNext())
						{
						String line = in.next();
						BedLine bed = this.bedCodec.decode(line);
						if(bed==null) continue;
						writer.write(line.getBytes());
						writer.write('\n');
						indexCreator.addFeature(bed, filePosition);
						filePosition = writer.getFilePointer();
						}
					}
				writer.flush();
				LOG.info("Creating index");
				final Index index = indexCreator.finalizeIndex(writer.getFilePointer());
				LOG.info("Writing index to "+tbi+ " using "+index.getClass());
				index.writeBasedOnFeaturePath(this.outputFile);
				}
			LOG.info("Done  N="+bedLineCount);
			} 
		catch (final Throwable e)
			{
			if(Files.exists(this.outputFile) && Files.isRegularFile(this.outputFile))
				{
				LOG.warning("Deleting "+this.outputFile);
				Files.delete(this.outputFile);
				if(Files.exists(tbi) && Files.isRegularFile(tbi)) {
					Files.delete(tbi);
					}
				}
			throw new IOException(e);
			}
		
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		LineIterator r = null;
		if(this.outputFile==null)
			{
			LOG.error(getProgramName()+" output file undefined.");
			return -1;
			}
		if(!this.outputFile.getFileName().toString().endsWith(".bed.gz"))
			{
			LOG.error("output file should en with .bed.gz but got "+outputFile);
			return -1;
			}
		try {
			if(args.isEmpty())
				{
				r = IOUtils.openStdinForLineIterator();
				}
			else if(args.size()==1)
				{
				r= IOUtils.openURIForLineIterator(args.get(0));
				}
			else
				{
				LOG.error("Illegal number of arguments.");
				return -1;
				}
			run(r);
			return 0;
		} catch (Exception e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			}
		}
	
	
	
		
		
	public static void main(String[] args)
		{
		new BedIndexTabix().instanceMainWithExit(args);
		}
	}
