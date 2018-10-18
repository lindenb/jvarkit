/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamFilterParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;


/**
BEGIN_DOC

## Input

input is a set of bam file or a file with the '*.list' suffix containing the path to the bam files.

## Example

TODO


END_DOC
 */
@Program(name="wescnvtview",
description="SVG visualization of bam DEPTH for multiple regions in a terminal",
keywords={"bam","alignment","graphics","visualization","svg","cnv"}
)
public class WesCnvTView  extends Launcher {
	private static final Logger LOG = Logger.build(WesCnvTView.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-B","--bed","-b","--capture"},description=
			"BED Capture. BED file containing the Regions to be observed.")
	private File bedFile = null;
	@Parameter(names={"-r","-rgn","--region"},description="Interval regions: 'CHR:START-END'. multiple separated with spaces or semicolon")
	private String bedRegions = null;
	@Parameter(names={"-w","--width","--cols","-C"},description="Terminal width")
	private int terminalWidth = 80 ;
	@Parameter(names={"-H","--height","--rows"},description="Terminal width")
	private int sampleHeight = 10 ;
	@Parameter(names={"-cap","--cap"},description="Cap coverage to this value. Negative=don't set any limit")
	private int capMaxDepth = -1 ;
	@Parameter(names={"--filter"},description=SamFilterParser.FILTER_DESCRIPTION,converter=SamFilterParser.StringConverter.class)
	private SamRecordFilter samRecordFilter = SamFilterParser.ACCEPT_ALL;
	@Parameter(names={"-p","-percentile","--percentile"},description="How to compute the percentil of a region")
	private Percentile.Type percentile = Percentile.Type.AVERAGE;
	
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
	private DecimalFormat decimalFormater = new DecimalFormat("##.##");
	private DecimalFormat niceIntFormat = new DecimalFormat("###,###");
	private final int LEFT_MARGIN = 6;		 
	
		 
	private abstract class AbstractViewWriter
		implements Closeable
		{
		int count_intervals = 0;
		final PrintWriter pw;
		final Map<String,SampleInfo> sample2info = new TreeMap<>();
		AbstractViewWriter(final PrintWriter pw) {
			this.pw = pw;
			}
		
		char float2val(double f)
			{
			char array[]= {' ','\u2581','\u2582','\u2583','\u2584','\u2585','\u2586','\u2587'};
			int idx =(int)(array.length*f);
			if(idx<=0) return ' ';
			if(idx>=array.length) return '\u2588';
			return array[idx];
			}
		
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
			this.sample2info.clear();
			this.pw.println(">>> " +labelOf(i));
			}
		void endInterval(final Interval i) {
			this.pw.println("<<< " +labelOf(i));
			this.sample2info.clear();
			}
		@Override
		public void close() throws IOException {
			this.pw.flush();
			this.pw.close();
			}
		void dump(final Interval interval) {
			int maxDepth = (int)this.sample2info.values().stream().
					flatMapToDouble(SI->Arrays.stream(SI.pixel_coverage)).
					max().orElse(0.0);
			if(WesCnvTView.this.capMaxDepth>0 && maxDepth>WesCnvTView.this.capMaxDepth) {
				maxDepth =  WesCnvTView.this.capMaxDepth;
				}
			if(maxDepth<=0) maxDepth=1;
			int idx=0;
			for(final String sample:this.sample2info.keySet())
				{
				if(idx>0)
					{
					this.pw.println();
					}
				dump(this.sample2info.get(sample),maxDepth,interval);
				++idx;
				}
			
			}
		void dump(final SampleInfo si,int maxDepth,final Interval interval) {
			LOG.debug(maxDepth);
			final int areaWidth = terminalWidth-LEFT_MARGIN;
			String s = "> "+si.sample+" ";
			this.pw.print(s);
			for(int i=s.length();i< terminalWidth;i++) this.pw.print("=");
			this.pw.println();
			//print ruler
				{
				s = "Pos: ";
				while(s.length() < LEFT_MARGIN)
					{
					s+=" ";
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
				}
			this.pw.println();
			for(int y=0;y< WesCnvTView.this.sampleHeight;++y)
				{
				double dpy = maxDepth - (maxDepth/(double) WesCnvTView.this.sampleHeight )*y;
				s= String.format("%.2f",dpy);
				while(s.length() < (LEFT_MARGIN-3))
					{
					s+=" ";
					}
				s+=" | ";
				pw.print(s);
				
				for(int x=0;x< si.pixel_coverage.length;x++)
					{
					
					if(si.pixel_coverage[x] < dpy)
						{
						pw.print(float2val(Math.abs(si.pixel_coverage[x]-Math.floor(si.pixel_coverage[x]))));
						}
					else
						{
						pw.print(float2val(1.0));
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
		}
		

	
	private void run(final AbstractViewWriter w,final Interval interval) {
			w.beginInterval(interval);
			for(final BamInput baminput: this.bamInputs)
				{
				run(w,baminput,interval);
				}
			w.dump(interval);
			w.endInterval(interval);
			w.sample2info.clear();
			}
	
	private void run(
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
			final Interval interval = new Interval(newCtg,interval0.getStart(),interval0.getEnd());
			final int base_coverage[] = new int[interval.length()];
			Arrays.fill(base_coverage, 0);
				
			final SAMRecordIterator iter=baminput.samReader.queryOverlapping(
						interval.getContig(),
						interval.getStart(),
						interval.getEnd()
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
				final Percentile thePercentile = Percentile.median();
				final int pos0 = Math.min(base_coverage.length, (int)(((x+0)/(double)si.pixel_coverage.length)*base_coverage.length));
				final int pos1 = Math.min(base_coverage.length, (int)Math.ceil(((x+1)/(double)si.pixel_coverage.length)*base_coverage.length));
				if(pos0>=pos1) continue;
				si.pixel_coverage[x] = thePercentile.evaluate(base_coverage,pos0,(pos1-pos0));
			}
			
		w.sample2info.put(baminput.sample,si);
		}

	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.terminalWidth-LEFT_MARGIN<10)
			{
			LOG.error("terminal width is too small");
			return -1;
			}
		BufferedReader r = null;
		PrintWriter out = null;
		try
			{			
			SAMSequenceDictionary firstDict = null;
			
			
			
			final SamReaderFactory srf = SamReaderFactory.makeDefault().
					validationStringency(ValidationStringency.LENIENT)
					;
			
			for(final File bamFile:IOUtils.unrollFiles2018(args)) {
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
			AbstractViewWriter w = new DefaultTerminalWriter(out);
			
			boolean got_interval = false;
			if(this.bedFile!=null)
				{
				String line;
				final BedLineCodec bedCodec = new BedLineCodec();
				r = IOUtils.openFileForBufferedReading(this.bedFile);
				while((line=r.readLine())!=null)
					{
					if(BedLine.isBedHeader(line)) continue;
					final BedLine bed = bedCodec.decode(line);
					if(bed==null || bed.getStart()>bed.getEnd()) {
						LOG.warn("Ignoring "+line);
						continue;
						}
					got_interval = true;
					run(w,bed.toInterval());
					if(out.checkError()) break;
					}
				r.close();
				r= null;
				}
			if(!StringUtil.isBlank(this.bedRegions))
				{
				final IntervalParser intervalParser = new IntervalParser(firstDict);
				for(final String s: this.bedRegions.split("[ \t;]+"))
					{
					if(StringUtil.isBlank(s)) continue;
					final Interval bed = intervalParser.parse(s);
					if(bed==null) {
						LOG.error("Cannot parse interval "+s);
						return -1;
						}
					run(w,bed);
					got_interval = true;
					if(out.checkError()) break;
					}
				}
			w.close();
			out.close();
			out = null;
			if(!got_interval) {
				LOG.warn("No interval was provided");
				}	
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			CloserUtil.close(r);
			CloserUtil.close(this.bamInputs);
			}
		}

	
public static void main(final String[] args)
	{
	new WesCnvTView().instanceMainWithExit(args);
	}
}
