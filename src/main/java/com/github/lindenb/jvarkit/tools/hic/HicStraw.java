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
package com.github.lindenb.jvarkit.tools.hic;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.List;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.hic.HicReader;
import com.github.lindenb.jvarkit.hic.HicReaderFactory;
import com.github.lindenb.jvarkit.hic.Normalization;
import com.github.lindenb.jvarkit.hic.Unit;
import com.github.lindenb.jvarkit.io.CustomSeekableStreamFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.seekablestream.ISeekableStreamFactory;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;

/**
BEGIN_DOC

### Example


```

END_DOC
 */

@Program(name="hicstraw",
	description="Query a Hi-C file",
	keywords={"hic"},
	creationDate="20190613",
	modificationDate="20190613",
	generate_doc=false
	)
public class HicStraw  extends Launcher {
	private static final Logger LOG = Logger.build(HicFileInfo.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-i","--interval1"},description="Interval 1",required=true)
	private String interval1Str = null;
	@Parameter(names={"-j","--interval2"},description="Interval 2. Use '*' to map all the chromosomes.",required=true)
	private String interval2Str = null;
	@Parameter(names={"-u","--unit"},description="Unit")
	private Unit unit = Unit.BP;
	@Parameter(names={"-n","--normalization"},description="normalization")
	private Normalization norm = Normalization.NONE;
	@Parameter(names={"-b","--bin"},description="bin size")
	private int binSize = 1;
	@Parameter(names={"-min-distance"},description="min distance between two intervals on the same chromosome. Don' print the value if they're closer than this value")
	private Integer minCisDistance = null;
	@Parameter(names={"-min-value"},description="Don' print the value if it's lower than 'v'")
	private Float minValue = null;
	@Parameter(names={"-max-value"},description="Don' print the value if it's greater than 'v'")
	private Float maxValue = null;

	
	
	private class MyQueryCallBack implements HicReader.QueryCallBack {
		PrintWriter pw;
		boolean first = true;
		String source;
		@Override
		public void reportContact(
				String contig1,int start1,int end1,
				String contig2,int start2,int end2,
				final Normalization norm,
				final Unit unit,
				final int binsize, 
				final float value
				)
			{
			if(this.first) {
				pw.println("##source="+this.source);
				pw.println("##unit="+unit);
				pw.println("##normalisation="+norm);
				pw.println("##bin-size="+binsize);
				pw.println("#CHROM1\tSTART1\tEND1\tCHROM2\tSTART2\tEND2\tVALUE");
				this.first = false;
				}
			if(minValue!=null && value < minValue.floatValue()) return;
			if(maxValue!=null && value > maxValue.floatValue()) return;
			
			if(minCisDistance!=null && contig1.equals(contig2)) {
				final int distance;
				if(CoordMath.overlaps(start1, end1, start2, end2)) {
					distance = 0;
					}
				else if(end1 < start2) {
					distance = start2 - end1;
					}
				else
					{
					distance = start1 - end2;
					}
				if(distance < minCisDistance) return;
				}
			pw.print(contig1);
			pw.print("\t");
			pw.print(start1);
			pw.print("\t");
			pw.print(end1);
			pw.print("\t");
			pw.print(contig2);
			pw.print("\t");
			pw.print(start2);
			pw.print("\t");
			pw.print(end2);
			pw.print(contig1);
			pw.print("\t");
			pw.print(start1);
			pw.print("\t");
			pw.print(value);
			pw.println();
			}
		};
	
		
		
	@Override
	public int doWork(final List<String> args) {
		try
			{
			final ISeekableStreamFactory seekableStreamFactory = new CustomSeekableStreamFactory();
			final MyQueryCallBack callback = new MyQueryCallBack();
			
			callback.pw= super.openPathOrStdoutAsPrintWriter(outputFile);
			
			for(final String input :args) {
				callback.source = input;
				
				try(final HicReader hicReader = new HicReaderFactory().
							setSeekableStreamFactory(seekableStreamFactory).
							open(input)) { 
				
					final Function<String,Locatable > parseInterval = (S)->{
						final Optional<Locatable> loc = hicReader.parseInterval(S);
						if(!loc.isPresent()) {
							LOG.error("bad interval : \""+S+"\" available are "+ hicReader.
									getDictionary().getSequences().stream().
									map(SR->SR.getSequenceName()).collect(Collectors.joining(" ; ")));
							return null;
							}
						return loc.get();
						};
					
					if(!hicReader.getBasePairResolutions().contains(this.binSize)) {
						LOG.error("bad binSize : \""+this.binSize+"\" available are "+ hicReader.getBasePairResolutions().stream().map(S->String.valueOf(S)).collect(Collectors.joining(" ; ")));
						return -1;
						}
						
					final Locatable loc1 = parseInterval.apply(this.interval1Str);
					if(loc1==null) return -1;
					
					final List<Locatable> loc2list;
					if("*".equals(this.interval2Str)) {
						final Locatable loc2 = parseInterval.apply(this.interval2Str);
						if(loc2==null) return -1 ;
						loc2list = java.util.Collections.singletonList(loc2);
						}
					else
						{
						loc2list = hicReader.getDictionary().
								getSequences().stream().
								map(SR->new Interval(SR.getSequenceName(),1,SR.getSequenceLength())).
								collect(Collectors.toList());
						}
					
					for(final Locatable loc2:loc2list) {
						callback.first = true;
						hicReader.query(loc1, loc2,norm, this.binSize, this.unit,callback);
						}
					}
				}
			callback.pw.flush();
			callback.pw.close();
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	
	public static void main(final String[] args) {
		new HicStraw().instanceMainWithExit(args);
		}
}
