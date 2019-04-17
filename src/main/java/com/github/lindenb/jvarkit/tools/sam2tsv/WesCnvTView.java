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
package com.github.lindenb.jvarkit.tools.sam2tsv;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamRecordFilterFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;


/**
BEGIN_DOC

## Input

Input is a set of regions to observe. It can be a 

   * BED file
   * VCF File with SV
   * Some intervals 'contig:start-end'

Input can be read from stdin

## Example

```
find src/test/resources/ -type f -name "S*.bam" > bam.list

$  java -jar dist/wescnvtview.jar -l bam.list  -P "RF01:100-200" 

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
java -jar dist/wescnvtview.jar --bams bam.list -P -F BED jeter.txt   |\
	csplit -b  '%05d.txt' -f cnv. -n 5 -s -z - '/^>>>/' '{*}' 
```

`-s`: "do not print counts of output file sizes". 

`-z`: "remove empty output files". 

`-n`: "use specified number of digits instead of 2". 

`-b`: "use sprintf FORMAT instead of %02d". 

`-f`: "prefix". 

## Note to self: view in less/more

```
java -jar dist/wescnvtview.jar --bams bam.list -P -F BED jeter.txt   | less -r
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
@Program(name="wescnvtview",
description="SVG visualization of bam DEPTH for multiple regions in a terminal",
keywords={"bam","alignment","graphics","visualization","svg","cnv"},
modificationDate="20190417"
)
public class WesCnvTView  extends Launcher {
	private static final Logger LOG = Logger.build(WesCnvTView.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-l","-B","--bams"},description=
			"The Bam file(s) to be displayed. If there is only one file which ends with '.list' it is interpreted as a file containing a list of paths")
	private List<File> theBamFiles = new ArrayList<>();

	@Parameter(names={"-w","--width","--cols","-C"},description="Terminal width. Under linux good idea is to use the environment variable ${COLUMNS}")
	private int terminalWidth = 80 ;
	@Parameter(names={"-H","--height","--rows"},description="Terminal width per sample")
	private int sampleHeight = 10 ;
	@Parameter(names={"-cap","--cap"},description="Cap coverage to this value. Negative=don't set any limit")
	private int capMaxDepth = -1 ;
	@Parameter(names={"--filter"},description=SamRecordFilterFactory.FILTER_DESCRIPTION,converter=SamRecordFilterFactory.class,splitter=NoSplitter.class)
	private SamRecordFilter samRecordFilter = SamRecordFilterFactory.getDefault();
	@Parameter(names={"-p","-percentile","--percentile"},description="How to compute the percentil of a region")
	private Percentile.Type percentile = Percentile.Type.MEDIAN;
	@Parameter(names={"-x","--extend"},description="Extend intervals by factor 'x'")
	private double extend_interval_factor = 1.0;
	@Parameter(names={"-F","--format"},description="input format. INTERVALS is a string 'contig:start-end'.")
	private InputFormat inputFormat = InputFormat.INTERVALS;
	@Parameter(names={"-P","--plain"},description="Plain output (not color)")
	private boolean plain_flag = false;
	@Parameter(names={"-highlight","--highlight","--top"},description="Per default samples are sorted alphabetically."
			+ "The samples in this collection will be displayed on the 'top' to have an quick insight about the propositus.")
	private Set<String> highlight_sample_set = new HashSet<>();
	@Parameter(names={"-G","--genes"},description="A BED file containing some regions of interest that will be displayed")
	private File roiFile = null;

	
	private enum InputFormat {VCF,BED,INTERVALS}
	
	private enum AnsiColor {
    	BLACK (30),
    	RED (31),
    	GREEN (32),
    	YELLOW (33),
    	BLUE (34),
    	MAGENTA (35),
    	CYAN (36),
    	WHITE (37)
		;
		final int opcode;
    	AnsiColor(final int opcode) {
    		this.opcode=opcode;
    		}
    	
     String color() {
    		return ANSI_ESCAPE+this.opcode+"m";
    		}
    	}

	
	public static final String ANSI_ESCAPE = "\u001B[";
	public static final String ANSI_RESET = ANSI_ESCAPE+"0m";
	private final IntervalTreeMap<Interval> geneMap = new IntervalTreeMap<>();
	private ContigNameConverter geneMapContigNameConverter = ContigNameConverter.getIdentity();

	
	private class BamInput implements Closeable
		{

		File bamFile;
		SamReader samReader=null;
		SAMSequenceDictionary dict = null;
		String sample;
		ContigNameConverter contigNameConverter;
		@Override
		public void close() throws IOException {
			CloserUtil.close(samReader);
			}
		}
	private static class SampleInfo
		{
		final String sample;
		double pixel_coverage[];
		SampleInfo(final String sample)
			{
			this.sample = sample;
			}
		}
	
	
	private final List<BamInput> bamInputs = new ArrayList<>();
	private DecimalFormat niceIntFormat = new DecimalFormat("###,###");
	private final int LEFT_MARGIN = 10;		 
	
		 
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
		
		protected String labelOf(final Interval i)
			{
			return i.getContig()+":"+
					niceIntFormat.format(i.getStart()) + "-" + 
					niceIntFormat.format(i.getEnd())+
					"\t Length:"+
					niceIntFormat.format(i.length()) +
					"\t("+count_intervals+")";
			}
			
		void beginInterval(final Interval i) {
			++count_intervals;
			if(count_intervals>1) this.pw.println();
			this.sampleInfos.clear();
			this.pw.println(">>> " +labelOf(i));
			}
		
		void printGenes(final Interval interval0) {
			final Interval interval = extendInterval(interval0);
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
					final int beg = interval.getStart()+ (int)((((x+0)-LEFT_MARGIN)/(double)(areaWidth))*interval.length());
					final int end = interval.getStart()+ (int)((((x+1)-LEFT_MARGIN)/(double)(areaWidth))*interval.length());
					
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
		
		void endInterval(final Interval i) {
			this.pw.println("<<< " +labelOf(i));
			this.sampleInfos.clear();
			}
		@Override
		public void close() throws IOException {
			this.pw.flush();
			this.pw.close();
			}
		void dump(final Interval interval0) {
			double maxDepth = this.sampleInfos.stream().
					flatMapToDouble(SI->Arrays.stream(SI.pixel_coverage)).
					max().orElse(0.0);
			
			if(WesCnvTView.this.capMaxDepth>0 && maxDepth>WesCnvTView.this.capMaxDepth) {
				maxDepth =  WesCnvTView.this.capMaxDepth;
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
		
		
		
		void dump(final SampleInfo si,double maxDepth,final Interval interval0)
			{
			final Interval interval = extendInterval(interval0);
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
					int pos = interval.getStart()+ (int)(((x-LEFT_MARGIN)/(double)(areaWidth))*interval.length());
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
			
			final double depthPerPixel = (1.0/ WesCnvTView.this.sampleHeight)*maxDepth;
			for(int pixy=0;pixy< WesCnvTView.this.sampleHeight;++pixy)
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
				
				
				final Function<Integer,Integer> pixel2base = (x)->interval.getStart()+(int)(((x)/(double)si.pixel_coverage.length)*interval.length());
				 
				
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
			super.pw.print(color.color());
			}
		@Override void endColor() {
			super.pw.print(ANSI_RESET);
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

	
	private void runInterval(final AbstractViewWriter w,final Interval interval) {
			w.beginInterval(interval);
			w.printGenes(interval);
			for(final BamInput baminput: this.bamInputs)
				{
				runInterval(w,baminput,interval);
				}
			w.dump(interval);
			w.endInterval(interval);
			w.sampleInfos.clear();
			}
	
	protected final Interval extendInterval(final Interval rgn) {
		final int x = (int)(rgn.length()*WesCnvTView.this.extend_interval_factor);
		if(x<=0) return rgn;
		return new Interval(
				rgn.getContig(),
				Math.max(rgn.getStart()-x,1),
				rgn.getEnd()+x
				);
		}
	
	private void runInterval(
			final AbstractViewWriter w,
			final BamInput baminput,
			final Interval interval0
			) {
			
			final SampleInfo si = new SampleInfo(baminput.sample);
			si.pixel_coverage = new double[this.terminalWidth-LEFT_MARGIN];
			Arrays.fill(si.pixel_coverage,0.0);
			
			final String newCtg = baminput.contigNameConverter.apply(interval0.getContig());
			if(StringUtil.isBlank(newCtg))
				{
				LOG.warn( "Contig not found in "+interval0+" for "+baminput.bamFile);
				return ;
				}
			final SAMSequenceRecord ssr = baminput.dict.getSequence(newCtg);
			if(ssr==null)
				{
				LOG.warn( "Contig not found in "+interval0+" for "+baminput.bamFile);
				return ;
				}
			
			final Interval interval1 = new Interval(newCtg,interval0.getStart(),interval0.getEnd());
					
			final Interval interval = extendInterval(interval1);
			
			final int base_coverage[] = new int[interval.length()];
			Arrays.fill(base_coverage, 0);
			
			final SAMRecordIterator iter=baminput.samReader.queryOverlapping(
						interval.getContig(),
						interval.getStart(),
						Math.min(ssr.getSequenceLength(),interval.getEnd())
						);
			while(iter.hasNext()) {
				final SAMRecord rec = iter.next();
				if(rec.getReadUnmappedFlag()) continue;
				if(this.samRecordFilter.filterOut(rec)) continue;
				
				final Cigar cigar=rec.getCigar();
				if(cigar==null || cigar.isEmpty()) continue;
				int ref1= rec.getStart();
					
				for(final CigarElement ce:cigar) {
					final CigarOperator op = ce.getOperator();					
					if(op.consumesReferenceBases())
						{
						if(op.consumesReadBases()){
							for(int x=0;x< ce.getLength();++x){
								final int pos=ref1+x;
								if(pos< interval.getStart()) continue;
								if(pos> interval.getEnd()) break;
								base_coverage[pos-interval.getStart()]++;
								}
							}
						ref1+=ce.getLength();
						}
					}
				}
			iter.close();
				
				
			for(int x=0;x< si.pixel_coverage.length;x++) {
				int pos0 = (int)(((x+0)/(double)si.pixel_coverage.length)*base_coverage.length);
				pos0 = Math.min(pos0,base_coverage.length);
				int pos1 = (int)Math.ceil(((x+1)/(double)si.pixel_coverage.length)*base_coverage.length);
				pos1 = Math.min(pos1,base_coverage.length);
				if(pos0>=pos1) continue;
				si.pixel_coverage[x] = Percentile.of(this.percentile).evaluate(base_coverage,pos0,(pos1-pos0));
				}
			
		w.sampleInfos.add(si);
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
	
	private Collection<Interval> getROIGenes(final Interval interval0) {
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
		PrintWriter out = null;
		try
			{
			final List<String> inputs = IOUtils.unrollStrings2018(args);
			
			
			final List<File> bamUrls = new ArrayList<>();
			if(this.theBamFiles.size()==1 && this.theBamFiles.get(0).getName().endsWith(".list"))
				{
				bamUrls.addAll(IOUtil.slurpLines(this.theBamFiles.get(0)).
						stream().
						filter(L->!StringUtil.isBlank(L) ).
						filter(L->!L.startsWith("#")).
						map(L->new File(L)).
						collect(Collectors.toSet())
						);
				}
			else
				{
				bamUrls.addAll(bamUrls);
				}
			
			if(bamUrls.isEmpty())
				{
				LOG.error("No BAM file was specified");
				return -1;
				}
			
			if(this.roiFile!=null) {
				BufferedReader br = IOUtils.openFileForBufferedReading(this.roiFile);
				final BedLineCodec codec = new BedLineCodec();
				String line;
				while((line=br.readLine())!=null)
					{
					final BedLine bed = codec.decode(line);
					if(bed==null) continue;
					final Interval interval = new Interval(
							bed.getContig(),
							bed.getStart(),
							bed.getEnd(),
							bed.getOrDefault(4, "+").equals("-"),
							bed.getOrDefault(3, "ROI."+(this.geneMap.size()+1))
							);
					this.geneMap.put(interval,interval);
					}
				br.close();
				this.geneMapContigNameConverter = ContigNameConverter.fromIntervalTreeMap(this.geneMap);
				}
			
			SAMSequenceDictionary firstDict = null;
			
			
			
			final SamReaderFactory srf = SamReaderFactory.makeDefault().
					validationStringency(ValidationStringency.LENIENT)
					;
			
			for(final File bamFile:bamUrls) {
				final BamInput bi = new BamInput();
				bi.bamFile = bamFile;
				bi.samReader = srf.open(bamFile);
				bi.dict = bi.samReader.getFileHeader().getSequenceDictionary();
				if(firstDict==null) firstDict = bi.dict;
				bi.contigNameConverter = ContigNameConverter.fromOneDictionary(bi.dict);
				JvarkitException.BamHasIndex.verify(bi.samReader);
				
				bi.sample =
					bi.samReader.getFileHeader().getReadGroups().stream().
					map(V->V.getSample()).
					filter(S->!StringUtil.isBlank(S)).
					findFirst().orElse(bamFile.getName());
					
				this.bamInputs.add(bi);
				}
			if(this.bamInputs.isEmpty()) {
				LOG.error("no bam input");
				return -1;
				}
			
			out = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			final AbstractViewWriter w = plain_flag?
					new PlainTerminalWriter(out):
					new DefaultTerminalWriter(out);
			
			switch(this.inputFormat)
				{
				case VCF:
					{
					final Predicate<VariantContext> acceptVariant = V->V.hasAttribute(VCFConstants.SVTYPE) && (V.hasAttribute("SVLEN") || V.hasAttribute("SVMETHOD")/* DELLY2 */) && V.getEnd()-V.getStart()>1;
					final Function<VariantContext,Interval> mapper = V->{
						int B = V.getStart();
						int E = V.getEnd();
						try {
							if(V.hasAttribute("CIPOS")) {
								final int x= V.getAttributeAsIntList("CIPOS", 0).get(0);
								B = Math.max(B-x,1);
								}
							if(V.hasAttribute("CIEND")) {
								final int x= V.getAttributeAsIntList("CIEND", 0).get(1);
								E +=x;
								}
							if(E<B) {
								final int tmp = B;
								B = E;
								E  = tmp;
								}
							}
						catch(final Throwable err) {
							
							}
						
						return new Interval(V.getContig(),B,E);
						};
					if(inputs.isEmpty())
						{
						final VCFIterator vcfin = super.openVCFIterator(null);
						while(vcfin.hasNext())
							{
							final VariantContext ctx = vcfin.next();
							if(!acceptVariant.test(ctx)) continue;
							runInterval(w,mapper.apply(ctx));
							if(out.checkError()) break;
							}
						vcfin.close();
						}
					else
						{
						for(final String vcfFile:inputs)
							{
							final VCFFileReader fr = new VCFFileReader(new File(vcfFile), false);
							fr.iterator().stream().
								filter(acceptVariant).
								map(mapper).
								forEach(I->runInterval(w,I));
							fr.close();
							if(out.checkError()) break;
							}
						}
					break;
					}
				case BED:
					{
					final BedLineCodec bedCodec = new BedLineCodec();
					final Consumer<BufferedReader> consummer = R->R.lines().
								filter(L->!StringUtil.isBlank(L)).
								map(L->bedCodec.decode(L)).
								filter(bed->!(bed==null || bed.getStart()>bed.getEnd())).
								map(B->B.toInterval()).
								forEach(I->{runInterval(w,I);})
								;					
					if(inputs.isEmpty())
						{
						final BufferedReader br = super.openBufferedReader(null);
						consummer.accept(br);
						br.close();
						}
					else
						{
						for(final String bedFile:inputs)
							{
							final BufferedReader br = IOUtils.openURIForBufferedReading(bedFile);
							consummer.accept(br);
							br.close();
							}
						}
					break;
					}
				case INTERVALS:
					{
					final IntervalParser parser=new IntervalParser();
					final Consumer<Stream<String>> consummer = SL->SL.filter(L->!StringUtil.isBlank(L)).
						map(L->parser.parse(L)).
						filter(I->I!=null && I.length()>1).
						forEach(I->{runInterval(w,I);})
						;
					if(inputs.isEmpty())
						{
						final BufferedReader br = super.openBufferedReader(null);
						consummer.accept(br.lines());
						br.close();
						}
					else
						{
						consummer.accept(inputs.stream());
						}
					break;
					}
				default: LOG.error("Invalid input format "+this.inputFormat); return -1;
				}
			
			w.close();
			out.close();
			out = null;
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			CloserUtil.close(this.bamInputs);
			}
		}

	
public static void main(final String[] args)
	{
	new WesCnvTView().instanceMainWithExit(args);
	}
}
