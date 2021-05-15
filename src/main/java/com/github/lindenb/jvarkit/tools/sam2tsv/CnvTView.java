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
package com.github.lindenb.jvarkit.tools.sam2tsv;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.ansi.AnsiUtils;
import com.github.lindenb.jvarkit.ansi.AnsiUtils.AnsiColor;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.CoverageFactory;
import com.github.lindenb.jvarkit.samtools.util.IntervalExtender;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;


/**
BEGIN_DOC

## Input

Input is a set of indexed cram/bam files or a file with the '.list' suffix containing the path to the bams

## Example

```
find src/test/resources/ -type f -name "S*.bam" > bam.list

$  java -jar dist/cnvtview.jar  -r "RF01:100-200" bam.list 

>>> RF01:100-200	 Length:101	(1)
> S1 ===========================================================================
     Pos| 1 9 18 31 44 56 69 82 95 108 125 142 160 177 194 211 228 246 263 280  
   9.00 |                                                                       
   8.10 |                                                                       
   7.20 |                %                                                      
   6.30 |           %%%%%%%%                                                    
   5.40 |           %%%%%%%%%      #           # #                              
   4.50 |       %%%%%%%%%%%%%   % ##  ####   #########   %%%%%                  
   3.60 |    %%%%%%%%%%%%%%%%%%%%########################%%%%%% %%%%   %%%     %
   2.70 |    %%%%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
   1.80 |   %%%%%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
   0.90 | %%%%%%%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%

> S2 ===========================================================================
     Pos| 1 9 18 31 44 56 69 82 95 108 125 142 160 177 194 211 228 246 263 280  
   9.00 |                                                           %%%  %%     
   8.10 |                                                           %%% %%%     
   7.20 |                                                   %%  %%%%%%%%%%%     
   6.30 |                     %                             %%%%%%%%%%%%%%%%%   
   5.40 |              %%%%%%%%   #       ###           #  %%%%%%%%%%%%%%%%%%   
   4.50 |          %%%%%%%%%%%%%%##  ##  #########    ###%%%%%%%%%%%%%%%%%%%%%%%
   3.60 |          %%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
   2.70 |      %%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
   1.80 |     %%%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
   0.90 |    %%%%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%

> S3 ===========================================================================
     Pos| 1 9 18 31 44 56 69 82 95 108 125 142 160 177 194 211 228 246 263 280  
   9.00 |                                                           %%%  %%     
   8.10 |                                                           %%% %%%     
   7.20 |                                                   %%  %%%%%%%%%%%     
   6.30 |                     %                             %%%%%%%%%%%%%%%%%   
   5.40 |              %%%%%%%%   #       ###           #  %%%%%%%%%%%%%%%%%%   
   4.50 |          %%%%%%%%%%%%%%##  ##  #########    ###%%%%%%%%%%%%%%%%%%%%%%%
   3.60 |          %%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
   2.70 |      %%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
   1.80 |     %%%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
   0.90 |    %%%%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%

> S4 ===========================================================================
     Pos| 1 9 18 31 44 56 69 82 95 108 125 142 160 177 194 211 228 246 263 280  
   9.00 |                                                                       
   8.10 |                                                                       
   7.20 |                                                                       
   6.30 |                                    ##                                 
   5.40 |                            ##########                                 
   4.50 |                      % ################        %  %%                  
   3.60 |                      %%################# ######%%%%%%%%%%             
   2.70 |                %%%%%%%%########################%%%%%%%%%%%%%          
   1.80 |        %%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%  %    
   0.90 |       %%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%

> S5 ===========================================================================
     Pos| 1 9 18 31 44 56 69 82 95 108 125 142 160 177 194 211 228 246 263 280  
   9.00 |                                                                       
   8.10 |                                                                       
   7.20 |                                                             %%       %
   6.30 |                                                           %%%%       %
   5.40 |                                                           %%%%%     %%
   4.50 |                       %#######                 %%%%%     %%%%%% %%%%%%
   3.60 |                   % %%%########                %%%%%%%%%%%%%%%%%%%%%%%
   2.70 |                %%%%%%%%########### #          #%%%%%%%%%%%%%%%%%%%%%%%
   1.80 |               %%%%%%%%%###############        #%%%%%%%%%%%%%%%%%%%%%%%
   0.90 |    %%%%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
<<< RF01:100-200	 Length:101	(1)


```

## Note to self: Splitting the output:

```
java -jar dist/cnvtview.jar --bams bam.list -P -F BED jeter.txt   |\
	csplit -b  '%05d.txt' -f cnv. -n 5 -s -z - '/^>>>/' '{*}' 
```

`-s`: "do not print counts of output file sizes". 

`-z`: "remove empty output files". 

`-n`: "use specified number of digits instead of 2". 

`-b`: "use sprintf FORMAT instead of %02d". 

`-f`: "prefix". 

## Note to self: view in less/more

```
java -jar dist/cnvtview.jar bam.list -r jeter.txt   | less -r
```

## Screenshot

https://twitter.com/yokofakun/status/1053185975923470337

![https://pbs.twimg.com/media/Dp2rDfsWoAEQEAI.jpg](https://pbs.twimg.com/media/Dp2rDfsWoAEQEAI.jpg)

https://twitter.com/yokofakun/status/1053204927202369536

![https://pbs.twimg.com/media/Dp28R1VWwAA7frV.jpg](https://pbs.twimg.com/media/Dp28R1VWwAA7frV.jpg)

https://twitter.com/yokofakun/status/1057627022665502721

![https://pbs.twimg.com/media/Dq1x60NVAAUhGPG.jpg](https://pbs.twimg.com/media/Dq1x60NVAAUhGPG.jpg)

END_DOC
 */
@Program(name="cnvtview",
description="Text visualization of bam DEPTH for multiple regions in a terminal",
keywords={"bam","alignment","graphics","visualization","cnv","ascii","text"},
modificationDate="20210412",
creationDate="20181018"
)
public class CnvTView  extends Launcher {
	private static final Logger LOG = Logger.build(CnvTView.class).make();
	private enum Format {plain,ansi};

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path reference = null;
	@Parameter(names={"-r","--regions","--interval"},description=IntervalListProvider.OPT_DESC,converter=IntervalListProvider.StringConverter.class,splitter=NoSplitter.class,required=true)
	private IntervalListProvider intervalListProvider = IntervalListProvider.empty();
	@Parameter(names={"-x","--extend"},description=IntervalExtender.OPT_DESC,converter=IntervalExtender.StringConverter.class,splitter=NoSplitter.class,required=false)
	private IntervalExtender intervalExtender = IntervalExtender.of("150%");
	@Parameter(names={"--mapq"},description="Min mapping quality")
	private int mappingQuality = 1;
	@Parameter(names={"-p","-percentile","--percentile"},description="How to compute the percentil of a region")
	private CoverageFactory.ScaleType percentile = CoverageFactory.ScaleType.MEDIAN;
	@Parameter(names={"--stddev"},description="Sort output on standard deviation")
	private boolean sort_on_stdev = false;

	
	@Parameter(names={"-w","--width","--cols","-C"},description="Terminal width. Under linux good idea is to use the environment variable ${COLUMNS}")
	private int terminalWidth = 80 ;
	@Parameter(names={"-H","--height","--rows"},description="Terminal width per sample")
	private int sampleHeight = 10 ;
	@Parameter(names={"-cap","--cap"},description="Cap coverage to this value. Negative=don't set any limit")
	private int capMaxDepth = -1 ;

	@Parameter(names={"--format"},description="output format")
	private Format output_format = Format.plain;
	@Parameter(names={"-highlight","--highlight","--top"},description="Per default samples are sorted alphabetically."
			+ "The samples in this collection will be displayed on the 'top' to have an quick insight about the propositus.")
	private Set<String> highlight_sample_set = new HashSet<>();
	@Parameter(names={"-G","--genes"},description="A BED file containing some regions of interest that will be displayed")
	private Path roiFile = null;
	@Parameter(names={"--flush"},description="do not wait for all bam to be scanned, do not sort, do not normalize on the depth of all bams, print the figure as soon as possible.")
	private boolean flushNow = false;


	
	private final IntervalTreeMap<Interval> geneMap = new IntervalTreeMap<>();
	private SAMSequenceDictionary dictionary = null;
	private ContigNameConverter geneMapContigNameConverter = ContigNameConverter.getIdentity();

	
	private static class SampleInfo
		{
		final String sample;
		double pixel_coverage[];
		SampleInfo(final String sample)
			{
			this.sample = sample;
			}
		double getStdDev() {
			if(pixel_coverage==null || pixel_coverage.length==0) return 0.0;
			final double avg = Arrays.stream(pixel_coverage).average().orElse(0.0);
			double n=0;
			for(double depth: pixel_coverage) {
				n+= Math.abs(avg  - depth);
				}
			return n/pixel_coverage.length;
			}
		}
	
	
	private final List<Path> bamInputs = new ArrayList<>();
	private final int LEFT_MARGIN = 10;		 
	private SamReaderFactory samReaderFactory=null;
		 
	private abstract class AbstractViewWriter
		implements Closeable
		{
		int count_intervals = 0;
		final PrintWriter pw;
		final List<SampleInfo> sampleInfos = new ArrayList<>();
		AbstractViewWriter(final PrintWriter pw) {
			this.pw = pw;
			}
		
		abstract void beginColor(AnsiColor color);
		abstract void endColor();
		abstract char float2val(double f);
		abstract void printBarSymbol(char c, boolean overlap_user_interval);
		abstract char genearrow(boolean negativestrand);
		
		protected String labelOf(final Locatable i)
			{
			return i.getContig()+":"+
					StringUtils.niceInt(i.getStart()) + "-" + 
					StringUtils.niceInt(i.getEnd())+
					"\t Length:"+
					StringUtils.niceInt(i.getLengthOnReference()) +
					"\t("+count_intervals+")";
			}
			
		void beginInterval(final Locatable i) {
			++count_intervals;
			if(count_intervals>1) this.pw.println();
			this.sampleInfos.clear();
			this.pw.println(">>> " +labelOf(i));
			}
		
		void printGenes(final Locatable interval0) {
			final Locatable interval = extendInterval(interval0);
			for(final Interval gene : getROIGenes(interval)) {
				String s = gene.getName();
				while(s.length() < (LEFT_MARGIN-3))
					{
					s =" "+s;
					}
				while(s.length() > (LEFT_MARGIN-3))
					{
					s = s.substring(1);
					}
				beginColor(AnsiColor.YELLOW);
				this.pw.print(s);
				endColor();
				this.pw.print(" | ");
				int x= LEFT_MARGIN;
				final int areaWidth = terminalWidth-LEFT_MARGIN;

				while(x<terminalWidth)
					{
					final int beg = interval.getStart()+ (int)((((x+0)-LEFT_MARGIN)/(double)(areaWidth))*interval.getLengthOnReference());
					final int end = interval.getStart()+ (int)((((x+1)-LEFT_MARGIN)/(double)(areaWidth))*interval.getLengthOnReference());
					
					if(CoordMath.overlaps(beg, end, gene.getStart(), gene.getEnd()))
						{
						this.pw.print(genearrow(gene.isNegativeStrand()));
						}
					else
						{
						this.pw.print(" ");
						}
					x++;
					}
				
				this.pw.println();
				}
		}
		
		void endInterval(final Locatable i) {
			this.pw.println("<<< " +labelOf(i));
			this.sampleInfos.clear();
			}
		@Override
		public void close() throws IOException {
			this.pw.flush();
			this.pw.close();
			}
		void dump(final Locatable interval0) {
			double maxDepth = this.sampleInfos.stream().
					flatMapToDouble(SI->Arrays.stream(SI.pixel_coverage)).
					max().orElse(0.0);
			
			if(CnvTView.this.capMaxDepth>0 && maxDepth > CnvTView.this.capMaxDepth) {
				maxDepth =  CnvTView.this.capMaxDepth;
				}
			if(maxDepth<=0.0) maxDepth=1.0;
			int idx=0;
			Collections.sort(this.sampleInfos,(A,B)->{
				final boolean topA = highlight_sample_set.contains(A.sample);
				final boolean topB = highlight_sample_set.contains(B.sample);
				if(topA!=topB) {
					if(topA) return -1;
					if(topB) return 1;
					}
				if (CnvTView.this.sort_on_stdev) {
					final double std1 = A.getStdDev();
					final double std2 = B.getStdDev();
					int i = Double.compare(std2, std1);//inverse sort
					if(i!=0) return i;
				}
				
				
				final int i=A.sample.compareTo(B.sample);
				return i;
				});
			for(final SampleInfo si:this.sampleInfos)
				{
				if(idx>0)
					{
					this.pw.println();
					}
				dump(si,maxDepth,interval0);
				++idx;
				}
			}
		
		
		
		void dump(final SampleInfo si,double maxDepth,final Locatable interval0)
			{
			final Locatable interval = extendInterval(interval0);
			final int areaWidth = terminalWidth-LEFT_MARGIN;
			String s = "> "+si.sample+" ";
			beginColor(AnsiColor.CYAN);
			this.pw.print(s);
			endColor();
			for(int i=s.length();i< terminalWidth;i++) this.pw.print("=");
			this.pw.println();
			//print ruler
				{
				beginColor(AnsiColor.GREEN);
				s = "Pos | ";
				while(s.length() < LEFT_MARGIN)
					{
					s= " "+s;
					}
				pw.print(s);	
				int x=s.length();
				while(x<terminalWidth)
					{
					int pos = interval.getStart()+ (int)(((x-LEFT_MARGIN)/(double)(areaWidth))*interval.getLengthOnReference());
					s = String.valueOf(pos);
					if(x+s.length()+1<terminalWidth)
						{
						s+=" ";
						}
					if(x+s.length()<terminalWidth)
						{
						pw.print(s);
						x+=s.length();
						}
					else
						{
						pw.print(" ");
						x++;
						}
					}
				endColor();
				this.pw.println();
				}
			
			//end print ruler
			
			final double depthPerPixel = (1.0/ CnvTView.this.sampleHeight)*maxDepth;
			for(int pixy=0;pixy< CnvTView.this.sampleHeight;++pixy)
				{
				final double dpy = maxDepth - pixy * depthPerPixel;
				s= String.format("%.2f",dpy);
				while(s.length() < (LEFT_MARGIN-3))
					{
					s =" "+s;
					}
				s+=" | ";
				if(dpy<20) beginColor(AnsiColor.YELLOW);
				pw.print(s);
				if(dpy<20) endColor();
				
				
				final Function<Integer,Integer> pixel2base = (x)->interval.getStart()+(int)(((x)/(double)si.pixel_coverage.length)*interval.getLengthOnReference());
				 
				
				for(int x=0;x< si.pixel_coverage.length;x++)
					{
					int chromStart = pixel2base.apply(x+0);
					int chromEnd = pixel2base.apply(x+1);
					final boolean overlap_extend = CoordMath.overlaps(chromStart, chromEnd, interval0.getStart(), interval0.getEnd());
					double depth = si.pixel_coverage[x];
					if(depth < dpy)
						{
						pw.print(' ');
						}
					else 
						{
						printBarSymbol(float2val(Math.abs(depth-dpy)/depthPerPixel), overlap_extend);
						}
					}
				pw.println();
				}
			
			}
		
		
		}
	private class DefaultTerminalWriter extends AbstractViewWriter
		{
		DefaultTerminalWriter(final PrintWriter pw) {
			super(pw);
			}
		@Override
		char genearrow(boolean negativestrand) {
			return negativestrand?'\u2190':'\u2192';
			}
		@Override 
		char float2val(double f)
			{
			char array[]= {' ','\u2581','\u2582','\u2583','\u2584','\u2585','\u2586','\u2587'};
			int idx =(int)(Math.ceil(array.length*f));
			if(idx<=0) return ' ';
			if(idx>=array.length) return '\u2588';
			return array[idx];
			}
		
		@Override void beginColor(final AnsiColor color) {
			super.pw.print(color.begin());
			}
		@Override void endColor() {
			super.pw.print(AnsiUtils.ANSI_RESET);
			}
		
		@Override
		void printBarSymbol(char c, boolean overlap_extend) {
			if(overlap_extend)
				{
				beginColor(AnsiColor.RED);
				}
			else
				{
				beginColor(AnsiColor.BLUE);
				}
			pw.print(c);
			endColor();
			}
		
		}
	private class PlainTerminalWriter extends AbstractViewWriter
		{
		PlainTerminalWriter(final PrintWriter pw) {
			super(pw);
			}
		@Override
		char genearrow(boolean negativestrand) {
			return negativestrand?'<':'>';
			}
		
		@Override 
		char float2val(double f)
			{
			return '#';
			}
		
		@Override void beginColor(final AnsiColor color) {
			}
		@Override void endColor() {
			}
		@Override
		void printBarSymbol(char c, boolean overlap_extend) {
			if(overlap_extend)
				{
				pw.print('#');
				}
			else
				{
				pw.print('%');
				}
			}
		}

	
	private void runInterval(final AbstractViewWriter w,final Locatable interval) {
			w.beginInterval(interval);
			w.printGenes(interval);
			int idx=0;
			for(final Path baminput: this.bamInputs)
				{
				try {
					final SampleInfo si = runInterval(w,baminput,interval);
					
					// flush now 
					if(this.flushNow) {
						if(idx>0) w.pw.println();
						double maxDepth = Arrays.stream(si.pixel_coverage).max().orElse(0.0);
						if(CnvTView.this.capMaxDepth>0 && maxDepth > CnvTView.this.capMaxDepth) {
							maxDepth =  CnvTView.this.capMaxDepth;
							}
						if(maxDepth<=0.0) maxDepth=1.0;
						w.dump(si,maxDepth,interval);
						++idx;
						}
					}
				catch(final IOException err) {
					throw new RuntimeIOException(err);
					}
				}
			if(!flushNow) w.dump(interval);
			w.endInterval(interval);
			w.sampleInfos.clear();
			}
	
	protected final Locatable extendInterval(final Locatable rgn) {
		return intervalExtender.apply(rgn);
		}
	
	private SampleInfo runInterval(
			final AbstractViewWriter w,
			final Path baminput,
			final Locatable interval0
			) throws IOException {
			
			
			final Locatable interval = extendInterval(interval0);
			
			try(SamReader samReader = samReaderFactory.open(baminput)) {
				final SAMFileHeader hdr = samReader.getFileHeader();
				SequenceUtil.assertSequenceDictionariesEqual(this.dictionary, hdr.getSequenceDictionary());
				JvarkitException.BamHasIndex.verify(samReader);
			
				final String sample = samReader.getFileHeader().getReadGroups().stream().
					map(V->V.getSample()).
					filter(S->!StringUtil.isBlank(S)).
					findFirst().orElse(IOUtils.getFilenameWithoutCommonSuffixes(baminput));
				
				final SampleInfo si = new SampleInfo(sample);
				si.pixel_coverage = new double[this.terminalWidth-LEFT_MARGIN];

				
				final CoverageFactory coverageFactory = new CoverageFactory().setMappingQuality(this.mappingQuality);
				final CoverageFactory.SimpleCoverage cov  = coverageFactory.getSimpleCoverage(samReader, interval, null);
				si.pixel_coverage = cov.scale(this.percentile, si.pixel_coverage.length);

				w.sampleInfos.add(si);
				return si;
				}
				
			}

	
	private static int getDefaultNumberOfColumns() {
		int ncols = -1;
		final String os= java.lang.System.getProperty("os.name","");
		if(os.toLowerCase().contains("linux")) {
			InputStream in = null;
			try
				{
				ProcessBuilder pb = new ProcessBuilder("bash","-c","tput cols").
							redirectErrorStream(true);
				final Process p = pb.start();
				in = p.getInputStream();
				final String colstr = IOUtils.copyToString(new InputStreamReader(in));
				in.close();
				in=null;
				p.waitFor();
				ncols = Integer.parseInt(colstr.trim());
				}
			catch(final Throwable err)
				{
				}
			}
		if(ncols<=0) ncols=80;
		return ncols;
		}
	
	private Collection<Interval> getROIGenes(final Locatable interval0) {
		if(this.geneMap.isEmpty()) return Collections.emptyList();
		final String contig = this.geneMapContigNameConverter.apply(interval0.getContig());
		if(StringUtil.isBlank(contig)) return Collections.emptyList();
		return this.geneMap.getOverlapping(new Interval(contig,interval0.getStart(),interval0.getEnd()));
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.terminalWidth<=0) {
			this.terminalWidth=getDefaultNumberOfColumns();
			}
		
		if(this.terminalWidth-LEFT_MARGIN<10)
			{
			LOG.error("terminal width is too small");
			return -1;
			}
		if(this.flushNow && !this.highlight_sample_set.isEmpty()) {
			LOG.warning("cannot highlight samples with option --flush.");
			}
		if(this.flushNow && this.sort_on_stdev) {
			LOG.warning("cannot sort on stddev option --flush.");
			}
		try
			{
			this.dictionary = SequenceDictionaryUtils.extractRequired(this.reference);
			final ContigNameConverter ctgConverter = ContigNameConverter.fromOneDictionary(this.dictionary);

			this.samReaderFactory = SamReaderFactory.makeDefault().
					referenceSequence(this.reference).
					validationStringency(ValidationStringency.LENIENT)
					;

			
			if(this.roiFile!=null) {
				try(BufferedReader br = IOUtils.openPathForBufferedReading(this.roiFile)) {
					final BedLineCodec codec = new BedLineCodec();
					String line;
					while((line=br.readLine())!=null)
						{
						final BedLine bed = codec.decode(line);
						if(bed==null) continue;
						final String ctg = ctgConverter.apply(bed.getContig());
						if(StringUtils.isBlank(ctg)) continue;
						final Interval interval = new Interval(
								ctg,
								bed.getStart(),
								bed.getEnd(),
								bed.getOrDefault(4, "+").equals("-"),
								bed.getOrDefault(3, "ROI."+(this.geneMap.size()+1))
								);
						this.geneMap.put(interval,interval);
						}
					}
				}
			
			this.bamInputs.addAll(IOUtils.unrollPaths(args));
			if(this.bamInputs.isEmpty()) {
				LOG.error("no bam input");
				return -1;
				}
			
			try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				final AbstractViewWriter w ;
				switch(this.output_format) {
						case plain: w=	new PlainTerminalWriter(out); break;
						case ansi: w = new DefaultTerminalWriter(out); break;
						default: throw new IllegalArgumentException();
					}
				this.intervalListProvider.
					dictionary(this.dictionary).
					stream().
					filter(R->!out.checkError()).
					forEach(R->runInterval(w,R));
				out.flush();
				}
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

	
public static void main(final String[] args)
	{
	new CnvTView().instanceMainWithExit(args);
	}
}
