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
package com.github.lindenb.jvarkit.tools.liftover;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.ucsc.LiftOverChain;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;

/**
BEGIN_DOC

## Example


```
cat in.bed | java -jar jvarkit.jar bedliftover --chain x.chain -R ref.fa
```


END_DOC
*/
@Program(
		name="bedliftover",
		description="LiftOver a BED file",
		creationDate="20140311",
		modificationDate="20240625",
		keywords={"bed","liftover"},
		jvarkit_amalgamion = true
		)
public class BedLiftOver extends Launcher
	{
	private static final Logger LOG = Logger.of(BedLiftOver.class);

	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;

	@Parameter(names={"-f","--chain"},description=LiftOverChain.OPT_DESC,required=true)
	private String liftOverFile = null;
	@Parameter(names={"-x","--failed"},description="  write bed failing the liftOver here. Optional.")
	private Path failedFile = null;
	@Parameter(names={"-m","--minmatch"},description="lift over min-match.")
	private double userMinMatch = LiftOver.DEFAULT_LIFTOVER_MINMATCH ;
	@Parameter(names={"-D","-R","-r","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null;
	@Parameter(names={"--chainvalid"},description="Ignore LiftOver chain validation")
	private boolean ignoreLiftOverValidation=false;
	@Parameter(names={"--original","--src"},description="Append original interval")
	private boolean appendOriginal=false;

	
	private void scan(final LiftOver liftOver,BufferedReader r,PrintWriter out,PrintWriter failed) throws IOException
		{
		String line;
		final BedLineCodec bedCodec=new BedLineCodec();
		while((line=r.readLine())!=null)
			{
			
			if(line.startsWith("#") || line.trim().isEmpty()) continue;
			final BedLine bedLine = bedCodec.decode(line);
			if(bedLine==null) continue;
			final Interval srcInterval = bedLine.toInterval();
			Interval dest=liftOver.liftOver(srcInterval);
			if(dest!=null)
				{
				out.print(dest.getContig());
				out.print('\t');
				out.print(dest.getStart()-1);
				out.print('\t');
				out.print(dest.getEnd());
				for(int i=3;i< bedLine.getColumnCount();++i) { 
					out.print('\t');
					out.print(bedLine.get(i));
					}
				if(this.appendOriginal) {
					out.print('\t');
					out.print(bedLine.getContig()+":"+bedLine.getStart()+"-"+bedLine.getEnd());
					}
				out.println();
				}
			else if(failed!=null)
				{
				failed.println(line);
				}			
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			final LiftOver liftOver= LiftOverChain.load(this.liftOverFile);
			liftOver.setLiftOverMinMatch(this.userMinMatch);

			
			if(!this.ignoreLiftOverValidation) {
				liftOver.validateToSequences(SequenceDictionaryUtils.extractRequired(faidx));
				}
			
			try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				try(PrintWriter failed=(this.failedFile!=null?super.openPathOrStdoutAsPrintWriter(this.failedFile):new PrintWriter(new NullOuputStream()))) {
					if(args.isEmpty())
						{
						try(BufferedReader r= openBufferedReader(null)) {
							scan(liftOver, r,out,failed);
							}
						}
					else
						{
						for(final String filename:args)
							{
							try(BufferedReader r=openBufferedReader(filename)) {
								scan(liftOver,r,out,failed);
								}
							}
						}
					
					failed.flush();
					}
				out.flush();
				}
			
			return 0;
			}
		catch(Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		}


	public static void main(final String[] args)
		{
		new BedLiftOver().instanceMainWithExit(args);
		}

	}
