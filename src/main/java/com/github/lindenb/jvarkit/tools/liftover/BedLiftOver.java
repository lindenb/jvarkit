/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLine;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.liftover.LiftOverChain;
import com.github.lindenb.jvarkit.liftover.LiftOverLoader;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

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
		modificationDate="20250526",
		keywords={"bed","liftover"},
		jvarkit_amalgamion = true
		)
public class BedLiftOver extends Launcher
	{
	private static final Logger LOG = Logger.of(BedLiftOver.class);

	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-c","--columns"},description="column indexes for chrom,start,end. Multiple chrom/start/end can be set by groups of 3 intergers: e.g '1,2,3,6,7,8,10,11,12' ")
	private String columnsStr = "1,2,3";
	@Parameter(names={"-f","--chain"},description=LiftOverChain.OPT_DESC,required=true)
	private String liftOverFile = null;
	@Parameter(names={"-x","--failed"},description="  write bed failing the liftOver here. Optional.")
	private Path failedFile = null;
	@Parameter(names={"-m","--minmatch"},description="lift over min-match.")
	private double userMinMatch = LiftOver.DEFAULT_LIFTOVER_MINMATCH ;
	@Parameter(names={"-R1"},description=LiftOverLoader.OPT_DICT_R1)
	private Path faidxR1 = null;
	@Parameter(names={"-D","-R2","-R","-r","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidxR2 = null;
	@Parameter(names={"--chainvalid","--disable-chain-validation"},description="Ignore LiftOver chain validation")
	private boolean ignoreLiftOverValidation=false;
	@Parameter(names={"--original","--src"},description="Append original interval as CHROM:START-END")
	private boolean appendOriginal=false;
	@Parameter(names={"-1"},description="coordinates are one-based (input is NOT bed)")
	private boolean one_based_index =false;

	
	private void scan(final LiftOver liftOver,
			final List<Integer> columns1,
			final BufferedReader r,PrintWriter out,PrintWriter failed) throws IOException
		{
		final int shift_bed = one_based_index?0:1;
		
		for(;;)
			{
			final String line = r.readLine();
			if(line==null) break;
			if( StringUtils.isBlank(line) || BedLine.isBedHeader(line)) continue;
			final String[] src = CharSplitter.TAB.split(line);
			final String[] tokens = Arrays.copyOf(src, src.length);
			boolean ok=true;
			final StringBuilder original = new StringBuilder();
			
			for(int k=0;k < columns1.size();k+=3) {
				final int contigIdx = columns1.get(k+0)-1;
				final int startIdx = columns1.get(k+1)-1;
				final int endIdx = columns1.get(k+2)-1;
				
				if(contigIdx<0 || contigIdx >= src.length ||
					startIdx< 0|| startIdx >= src.length ||
					endIdx<0 || endIdx >= src.length) {
					ok=false;
					LOG.warn("not enough columns in $"+(contigIdx+1)+" $"+(startIdx+1)+" $"+(endIdx+1)+":\""+String.join("(tab)", src)+"\"");
					break;
					}
				final String contig =  src[contigIdx];
				final int start1 = Integer.parseInt(src[startIdx]) + shift_bed;//bed to interval
				if(start1 < (one_based_index?1:0)) {
					LOG.warn("Bad start in :\""+String.join("(tab)", src)+"\"");
					break;
					}
				final int end1 = Integer.parseInt(src[endIdx]);
				
				final Interval srcInterval = new Interval(contig,start1,end1);
				final Interval dest=liftOver.liftOver(srcInterval);
				if(dest==null) {
					ok=false;
					break;
					}
				tokens[contigIdx]= dest.getContig();
				tokens[startIdx]= String.valueOf(dest.getStart()-shift_bed);
				tokens[endIdx]= String.valueOf(dest.getEnd());
				
				if(this.appendOriginal) {
					if(original.length()>0) original.append(",");
					original.append(contig+":"+start1+"-"+end1);
					}
				} // end loop over columns
			
			if(ok) {
				out.print(String.join("\t", tokens));
				if(this.appendOriginal) {
					out.print('\t');
					out.print(original.toString());
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
	public int doWork(final List<String> args0) {
		try
			{
			final List<String> args  = IOUtils.unrollStrings(args0);
			List<Integer> columns= new ArrayList<>();
			for(String s:CharSplitter.COMMA.split(this.columnsStr)) {
				try {
					final int i = Integer.parseInt(s);
					if(i<1) {
						LOG.error("bad column index in "+this.columnsStr);
						return -1;
						}
					columns.add(i);
					}
				catch(NumberFormatException err) {
					LOG.error(err);
					return -1;
				}
			}
			
			if(columns.isEmpty()) {
				LOG.error("empty columns");
				return -1;
				}
			if(columns.size()%3!=0) {
				LOG.error("expected a multiple of 3 columns in "+columns);
				return -1;
				}
			final LiftOver liftOver= new LiftOverLoader().
					setSourceMapper(this.faidxR1==null?null:ContigNameConverter.fromPathOrOneDictionary(this.faidxR1)).
					setTargetMapper(this.faidxR2==null?null:ContigNameConverter.fromPathOrOneDictionary(this.faidxR2)).
					load(this.liftOverFile);
			liftOver.setLiftOverMinMatch(this.userMinMatch);


			if(!this.ignoreLiftOverValidation) {
				liftOver.validateToSequences(SequenceDictionaryUtils.extractRequired(faidxR2));
				}

			try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {

				try(PrintWriter failed=(this.failedFile!=null?super.openPathOrStdoutAsPrintWriter(this.failedFile):NullOuputStream.newPrintWriter())) {
					if(args.isEmpty())
						{
						try(BufferedReader r= IOUtils.openStreamForBufferedReader(stdin())) {
							scan(liftOver,columns, r,out,failed);
							}
						}
					else
						{
						for(final String filename:args)
							{
							try(BufferedReader r= IOUtils.openURIForBufferedReading(filename)) {
								scan(liftOver,columns,r,out,failed);
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
