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

*/
package com.github.lindenb.jvarkit.tools.liftover;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.CloserUtil;


@Program(
		name="bedliftover",
		description="Lift-over a VCF file",
		keywords={"bed","liftover"}
		)
public class BedLiftOver extends Launcher
	{
	private static final Logger LOG = Logger.build(BedLiftOver.class).make();

	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-f","--chain"},description="LiftOver file.",required=true)
	private File liftOverFile = null;
	@Parameter(names={"-x","--failed"},description="  write bed failing the liftOver here. Optional.")
	private File failedFile = null;
	@Parameter(names={"-m","--minmatch"},description="lift over min-match.")
	private double userMinMatch = LiftOver.DEFAULT_LIFTOVER_MINMATCH ;
	@Parameter(names={"-D","-R","-r","--reference"},description="indexed REFerence file.",required=true)
	private File faidx = null;
	@Parameter(names={"--chainvalid"},description="Ignore LiftOver chain validation")
	private boolean ignoreLiftOverValidation=false;
	
	private LiftOver liftOver=null;

	
	private void scan(BufferedReader r,PrintWriter out,PrintWriter failed) throws IOException
		{
		String line;
		final BedLineCodec bedCodec=new BedLineCodec();
		while((line=r.readLine())!=null)
			{
			
			if(line.startsWith("#") || line.trim().isEmpty()) continue;
			final BedLine bedLine = bedCodec.decode(line);
			if(bedLine==null) continue;
			final Interval srcInterval = bedLine.toInterval();
			Interval dest=this.liftOver.liftOver(srcInterval);
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
				out.println();
				}
			else if(failed!=null)
				{
				failed.println(line);
				}			
			}
		}
	
	@Override
	public int doWork(List<String> args) {
		
		if(liftOverFile==null)
			{
			LOG.error("LiftOver file is undefined.");
			return -1;
			}
		this.liftOver=new LiftOver(liftOverFile);
		this.liftOver.setLiftOverMinMatch(this.userMinMatch);
		
		PrintWriter out=null;
		PrintWriter failed=null;
		try
			{
			if(!this.ignoreLiftOverValidation) {
				IndexedFastaSequenceFile ref=new IndexedFastaSequenceFile(faidx);
				this.liftOver.validateToSequences(ref.getSequenceDictionary());
				ref.close();
				}
			
			out = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			
			if(this.failedFile!=null)
				{
				failed= super.openFileOrStdoutAsPrintWriter(failedFile);
				}
			if(args.isEmpty())
				{
				BufferedReader r= openBufferedReader(null);
				scan(r,out,failed);
				CloserUtil.close(r);
				}
			else
				{
				for(final String filename:args)
					{
					BufferedReader r=openBufferedReader(filename);
					scan(r,out,failed);
					CloserUtil.close(r);
					}
				}
			out.flush();
			out.close();
			out=null;
			if(failed!=null) {
				failed.flush();
				failed.close();
				failed=null;
			}
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			CloserUtil.close(failed);
			}
		}


	public static void main(final String[] args)
		{
		new BedLiftOver().instanceMainWithExit(args);
		}

	}
