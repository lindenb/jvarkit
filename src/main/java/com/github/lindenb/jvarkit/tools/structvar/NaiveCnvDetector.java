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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.iterator.MergingIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;


/**
BEGIN_DOC

Input is either:

  * one fileof samtools depth output. All 'N' samples in one file.
  * 'N' files (samtools depth output AND/OR bigwig/bigbed [experimental not tested] ). One samples in per file. REF dictionary is required. List of file can be specified if input ends with '.list' 

bigbed and bigwig have not been tested; Bigbed shouldn't have overlapping regions...

## Example

```
samtools depth -r '1:1234-567' *.bam |\
	java -jar dist/naivecnvdetector.jar  > out.tsv
```


END_DOC

 */
@Program(name="naivecnvdetector",
	description="experimental CNV detection for multiple samples.",
	keywords= {"cnv","bam","sam","wig","bigwig","bigbed"}
	)
public class NaiveCnvDetector extends Launcher
	{
	private static final Logger LOG = Logger.build(NaiveCnvDetector.class).make();
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@Parameter(names={"-c","--config"},description="config file. Tab delimited. Sample-name(tab)mean-depth(tab)integer[affected=1,non-affected=0]. If this file is not specified , all samples are considered unaffected (discovery mode).")
	private File configFile=null;	
	/** size of a window */
	@Parameter(names={"-w"},description="window size")
	private int windowSize=1000;
	@Parameter(names={"-s"},description="window shift")
	private int windowShift=500;

	@Parameter(names={"--weirdDepth"},description="Treat normalized depth greater than this value as 'weird' and discard the sliding windows at this place.")
	private int weirdDepth=500;
	@Parameter(names={"-md","--min-dp"},description="At least one 'unaffected' sample must have a normalized-depth greater than this value.")
	private int min_unaffected_depth = 20;
	@Parameter(names={"--disable-consecutive"},description="Disable 'consecutive' positions criteria. Default: dump buffer of position if the current samtools-depth line is not the very next expected position.")
	private boolean disable_consecutive_bases =false;
	@Parameter(names={"-stddevu","--stddev-unaffected"},description="Maximum standard deviation of depth for unaffected samples. Ignored if negative or if not any affected samples is defined." )
	private double max_stdev_unaffected = 10.0;
	@Parameter(names={"-E","-del","--del","--deletion"},description="Deletion Treshold. Which fraction of the median depth is considered as aa deletion. Must be <1.0" )
	private double deletion_treshold = 0.5;
	@Parameter(names={"-U","-dup","--dup","--duplication"},description="Duplication Treshold. Which fraction of the median depth is considered as a duplication. Must be >1.0" )
	private double duplication_treshold = 1.5;
	@Parameter(names={"--no-both"},description="There cannot be a DEL and a DUP at the same place." )
	private boolean no_both = false;
	@Parameter(names={"-t"},description="DEL must be < median-depth-stdev and DUP must be > median-depth+stdev" )
	private boolean use_standard_depth = false;
	@Parameter(names={"-R","-reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private File dictRefFile =  null;



	
	private class SampleInfo
		{
		String name;
		int index;
		double meanDepth=0.0;
		double adjustDepth = 1.0;
		boolean affected=false;
		long sumDepth = 0L;
		long countDepth = 0L;
		
		boolean isAffected() { return affected;}
		boolean isUnaffected() { return !affected;}
		
		String getLabel() {
			return name+(count_affected_samples>0?(isAffected()?"*":""):"");
		}
		
		@Override
		public int hashCode() {
			return Integer.hashCode(this.index);
			}
		@Override
		public boolean equals(Object obj) {
			return SampleInfo.class.cast(obj).index==this.index;
			}
		@Override
		public String toString() {
			return name+"["+(index+1)+"]";
			}
		}
	
	private final List<SampleInfo> sampleList  = new ArrayList<>();
	private int count_affected_samples = 0;
	private int count_unaffected_samples = 0;
	private DecimalFormat decimalFormater = new DecimalFormat("##.##");
	private SAMSequenceDictionary dict = null;
	
	private class DepthLine
		{
		final String contig;
		final int pos;
		final double depths[];
		DepthLine(final String contig,final int pos,int num_samples){
			this.contig = contig;
			this.pos = pos;
			this.depths = new double[num_samples];
			Arrays.fill(this.depths, 0);
			}
		}
	private final int smooth_win = 10;
	private class DepthInterval
		{
		final String contig;
		final int start;
		final int end;
		final double depths[];
		DepthInterval(final List<DepthLine> lines){
			this.contig = lines.get(0).contig;
			this.start = lines.get(0).pos;
			this.end = lines.get(lines.size()-1).pos;
			final int n_samples= lines.get(0).depths.length;
			this.depths = new double[n_samples];
			Arrays.fill(this.depths, 0.0);
			double one_sample_depth[]=new double[lines.size()];
			double smooth_sample_depth[]=new double[lines.size()];
			for(int sample_index=0;sample_index< n_samples;sample_index++) {
				Arrays.fill(one_sample_depth,0);
				for(int j=0;j<lines.size();++j)
					{
					one_sample_depth[j]=lines.get(j).depths[sample_index];
					}
				//smooth values
				if(smooth_win>0) {
					Arrays.fill(smooth_sample_depth,0);
					for(int x=0;x< one_sample_depth.length;x++)
						{
						double t=0;
						int N=0;
						for(int y=-smooth_win;y<=smooth_win ;++y)
							{
							if(x+y<0) continue;
							if(x+y>=one_sample_depth.length) break;
							t+=one_sample_depth[x+y];
							N++;
							}
						smooth_sample_depth[x]=t/N;
						}
					System.arraycopy(smooth_sample_depth, 0, one_sample_depth, 0, smooth_sample_depth.length);
					}
				this.depths[sample_index] = new Mean().evaluate(one_sample_depth);
				}
			}
		
		}
	
	private final List<DepthLine> depthBuffer=new ArrayList<>();
	
	/** convert double to string */
	private String format(double v)
		{
		return this.decimalFormater.format(v);
		}	
	
	private void dump(final PrintWriter out) {
		if(depthBuffer.isEmpty()) return;
		if(depthBuffer.size()< this.windowSize/2) return;
		final DepthInterval rec = new DepthInterval(this.depthBuffer);
		
		// at last one unaffected must have depth >= this.min_unaffected_depth
		if(this.sampleList.stream().
			filter(S->S.isUnaffected()).
			noneMatch(S->rec.depths[S.index]>= this.min_unaffected_depth))
			{
			return;
			}
		
		//any sample having a weird depth
		if(Arrays.stream(rec.depths).anyMatch(V->V>=weirdDepth)) return;

		
		// standard deviation of unaffected
		final double stddev_unaffected = new StandardDeviation(true).
				evaluate(
				this.sampleList.stream().
				filter(S->S.isUnaffected()).
				mapToDouble(S->rec.depths[S.index]).
				toArray()
				);
		
		// there are some affected sample but stddev of unaffected is large
		if(this.count_affected_samples>0 && 
			this.max_stdev_unaffected>=0 && 
			stddev_unaffected>this.max_stdev_unaffected)
			{
			return;
			}
		
		
		
		//calc median depth of unaffected
		final double median_unaffected_depth = new Median().
				evaluate(
				this.sampleList.stream().
				filter(S->S.isUnaffected()).
				mapToDouble(S->rec.depths[S.index]).
				toArray()
				);
		
		if(median_unaffected_depth<=0 || Double.isNaN(median_unaffected_depth)) return;
		
		
		final Predicate<SampleInfo> deletionTest = (SI)->{
			final double dp = rec.depths[SI.index];
			if(use_standard_depth && dp> median_unaffected_depth-stddev_unaffected) return false;
			return dp <= median_unaffected_depth*this.deletion_treshold;
			};
			
		final Predicate<SampleInfo> duplicationTest = (SI)->{
			final double dp = rec.depths[SI.index];
			if(use_standard_depth && dp< median_unaffected_depth+stddev_unaffected) return false;
			return dp >= median_unaffected_depth*this.duplication_treshold;
			};	
		
		final List<SampleInfo> delSamples = this.sampleList.
				stream().
				filter(deletionTest).
				collect(Collectors.toList())
				;
		final List<SampleInfo> dupSamples = this.sampleList.
				stream().
				filter(duplicationTest).
				collect(Collectors.toList())
				;
		
		if(delSamples.isEmpty() && dupSamples.isEmpty()) return;
		
		// at least one sample affected must carry mutation
		if(this.count_affected_samples>0) {
			if(delSamples.stream().noneMatch(S->S.isAffected()) &&
				dupSamples.stream().noneMatch(S->S.isAffected()))
				{
				return;
				}
		}
		
		
		final Set<SampleInfo> noCnvSamples = new HashSet<>(this.sampleList);
		noCnvSamples.removeAll(delSamples);
		noCnvSamples.removeAll(dupSamples);
		
		
		//interval contains DEL *and* DUP
		if(this.no_both && 
			!delSamples.isEmpty() && 
			!dupSamples.isEmpty()
			)
			{
			return;
			}
		out.print(rec.contig);
		out.print("\t");
		out.print(rec.start-1);
		out.print("\t");
		out.print(rec.end);
		out.print("\t");
		out.print(format(median_unaffected_depth));
		out.print("\t");
		out.print(format(stddev_unaffected));
		
		for(int side=0;side<2;++side)
			{
			final List<SampleInfo> list = side==0?delSamples:dupSamples;
			out.print("\t");
			out.print(list.isEmpty()?".":(side==0?"DEL":"DUP"));
			out.print("\t");
			out.print(list.size());
			if(this.count_affected_samples>0) {
				out.print("\t");
				out.print(list.stream().filter(S->S.isAffected()).count());
				out.print("\t");
				out.print(list.stream().filter(S->S.isUnaffected()).count());
				}
			
			
			final DoubleStream dst = list.stream().
						mapToDouble(si->rec.depths[si.index]/median_unaffected_depth);
			final OptionalDouble best_ratio=(side==0?dst.min():dst.max());
				
			out.print("\t");
			out.print(best_ratio.isPresent()?format(best_ratio.getAsDouble()):".");
			out.print("\t");
			out.print(list.isEmpty()?".":list.stream().map(S->S.getLabel()).collect(Collectors.joining(";")));
			}
		
		
		
		if(this.count_affected_samples>0)
			{
			final ChiSquareTest chiSquareTest = new ChiSquareTest();
			
			final double p_value=chiSquareTest.chiSquare(new long[][]{
				{
				dupSamples.stream().filter(S->S.isAffected()).count() + delSamples.stream().filter(S->S.isAffected()).count(),
				dupSamples.stream().filter(S->S.isUnaffected()).count() + delSamples.stream().filter(S->S.isUnaffected()).count(),
				},{
				noCnvSamples.stream().filter(S->S.isAffected()).count(),
				noCnvSamples.stream().filter(S->S.isUnaffected()).count()
				}}
				);
			out.print("\t");
			out.print(format(p_value));
			out.print("\t");
			out.print(p_value<0.05?"*":".");
			}
		
		for(final SampleInfo si: this.sampleList)
			{
			out.print("\t");
			out.print(format(rec.depths[si.index])+" x"+format(rec.depths[si.index]/median_unaffected_depth));
			}
		
		out.println();
		out.flush();
		}
	
	
	private void printHeader(final PrintWriter out) {
		out.print("#chrom");
		out.print("\t");
		out.print("start");
		out.print("\t");
		out.print("end");
		out.print("\t");
		out.print("median.unaffected.depth");
		out.print("\t");
		out.print("median.unaffected.depth.stddev");
		for(int side=0;side<2;++side)
			{
			final String prefix=side==0?"DEL":"DUP";
			out.print("\t");
			out.print(prefix);
			out.print("\t");
			out.print(prefix+".count");
			if(this.count_affected_samples>0) {
				out.print("\t");
				out.print(prefix+".count.affected");
				out.print("\t");
				out.print(prefix+".count.unaffected");
				}
			out.print("\t");
			out.print(prefix+".best.ratio");
			out.print("\t");
			out.print(prefix+".samples");
			}
		
		
		if(this.count_affected_samples>0) {
			out.print("\t");
			out.print("chi2");
			out.print("\t");
			out.print("chi2.signifiant");
			}
		for(final SampleInfo si: this.sampleList)
			{
			out.print("\t");
			out.print(si.getLabel());
			}
		
		out.println();
		}

	
	@SuppressWarnings("resource")
	@Override
	public int doWork(final List<String> args) {		
	
		if(this.windowSize<10) {
			LOG.error("low window size");
			return -1;
		}
		if(this.windowShift<10) {
			LOG.error("low window shift");
			return -1;
		}
		if(this.deletion_treshold>=1.0)
			{
			LOG.error("bad deletion treshold . Must be <1.0 but got "+this.deletion_treshold);
			return -1;
			}
		if(this.duplication_treshold<=1.0)
			{
			LOG.error("bad dup treshold . Must be >1.0 but got "+this.duplication_treshold);
			return -1;
			}
		
		
		
		PrintWriter out = null;
		try
			{
			if(this.dictRefFile!=null) 
				{
				this.dict = SAMSequenceDictionaryExtractor.extractDictionary(this.dictRefFile);
				}
			else
				{
				this.dict=null;
				}
			
			
			out =  super.openFileOrStdoutAsPrintWriter(this.outputFile);

			
			if(this.configFile!=null)
				{
				this.sampleList.addAll( IOUtil.slurpLines(this.configFile).stream().
						filter(S->!(StringUtil.isBlank(S) || S.startsWith("#"))).
						map(S->JvarkitException.TokenErrors.atLeast(3,CharSplitter.TAB.split(S))).
						map(T->{
							final SampleInfo si=new SampleInfo();
							si.name = T[0];
							si.meanDepth = Double.parseDouble(T[1]);
							if(si.meanDepth<0) throw new IllegalArgumentException("negative mean depth");
							si.affected=false;
							if(T[2].equalsIgnoreCase("true") || T[2].equals("1")) {
								si.affected=true;
								}
							return si;
							}).collect(Collectors.toList())
						);
				for(int i=0;i< sampleList.size();i++)
					{
					this.sampleList.get(i).index=i;
					}
				final double max_depth= this.sampleList.stream().mapToDouble(S->S.meanDepth).max().getAsDouble();
				for(int i=0;i< this.sampleList.size();i++)
					{
					this.sampleList.get(i).adjustDepth= max_depth / this.sampleList.get(i).meanDepth ;
					LOG.debug("Max-depth "+max_depth+" "+this.sampleList.get(i).meanDepth+" "+this.sampleList.get(i).adjustDepth);
					}
				this.count_affected_samples = (int)this.sampleList.stream().filter(S->S.affected).count();
				this.count_unaffected_samples = this.sampleList.size()-this.count_affected_samples;
				
				if(this.count_unaffected_samples==0)
					{
					LOG.error("No unaffected sample in config file "+this.configFile);
					return -1;
					}
				LOG.info("affected "+this.count_affected_samples+" ; unaffected:"+this.count_unaffected_samples);
				printHeader(out);
				}
			
			final Iterator<DepthLine> dpIter;
			if(args.isEmpty() || (args.size()==1 && !args.get(0).endsWith(".list"))) {
				dpIter = new MultipleSampleDepthIterators(super.openBufferedReader(oneFileOrNull(args)));
				}
			else
				{
				if(this.dict==null) {
					LOG.error("Ref dictionary is required . args: "+String.join(",", args)+".");
					return -1;
					}
				final List<File> filenames = IOUtils.unrollFiles2018(args);	
					if(filenames.isEmpty()) {
					LOG.error("no input");
					return -1;
					}
				dpIter = new CombineOneSampleDepthIterators(this.dict,filenames);
				}
			
			while(dpIter.hasNext()) {
				final DepthLine depthLine = dpIter.next();
				if(this.sampleList.isEmpty())
					{
					LOG.info("building default sample list");
					this.count_affected_samples = 0;
					this.count_unaffected_samples = depthLine.depths.length;
					for(int i=0;i< this.count_unaffected_samples;i++) {
						final SampleInfo si = new SampleInfo();
						si.index=i;
						si.affected=false;
						si.meanDepth = 20;
						si.adjustDepth = 1.0;
						si.name = String.format("S%03d",(i+1));
						this.sampleList.add(si);
						}
					printHeader(out);
					}
				else if(this.sampleList.size()!= depthLine.depths.length)
					{
					LOG.error(
							"expected at least "+( this.sampleList.size()+3)+" words"+depthLine.depths.length);
					return -1;
					}
					
				for(int x=0;x<depthLine.depths.length;++x)
					{
					final SampleInfo si=this.sampleList.get(x);					
					si.sumDepth += depthLine.depths[x];
					si.countDepth++;
					depthLine.depths[x] *= si.adjustDepth;
					}
				if(!this.depthBuffer.isEmpty())
					{
					final DepthLine last = this.depthBuffer.get(this.depthBuffer.size()-1);
					if(!last.contig.equals(depthLine.contig))
						{
						dump(out);
						this.depthBuffer.clear();
						}
					else if(last.pos+1!=depthLine.pos)
						{
						dump(out);
						this.depthBuffer.clear();
						}
					else if(!this.disable_consecutive_bases && last.pos>=depthLine.pos)
						{
						dump(out);
						this.depthBuffer.clear();
						}
					}
				this.depthBuffer.add(depthLine);
				if(this.depthBuffer.size()==this.windowSize)
					{
					dump(out);
					for(int x=0;x<this.windowShift && !this.depthBuffer.isEmpty();++x)
						{
						this.depthBuffer.remove(0);
						}
					}
				}

			out.flush();
			out.close();
			CloserUtil.close(dpIter);
			
			for(final SampleInfo ci:this.sampleList) {
				if(ci.countDepth<=0) continue;
				LOG.info(ci.name+"\t"+(ci.sumDepth/(double)ci.countDepth)+"\t"+(ci.isAffected()?1:0));
				}
			
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			}
		}
	
	private static class OneSampleDepth
			{
			final int sample_index ;
			final int tid;
			final String contig;
			final int pos1;
			final int rawdp;
			OneSampleDepth(final int sample_index,final SAMSequenceDictionary dict,final String line) {
				this.sample_index = sample_index;
				final String tokens[] = CharSplitter.TAB.split(line);
				if(tokens.length!=3) {
					throw new JvarkitException.TokenErrors("expected  3 words",tokens);
					}
				this.contig = tokens[0];
				this.tid = dict.getSequenceIndex(contig);
				if(this.tid==-1) throw new JvarkitException.ContigNotFoundInDictionary(this.contig, dict);
				this.pos1 = Integer.parseInt(tokens[1]);
				this.rawdp = Integer.parseInt(tokens[2]);
				}
			OneSampleDepth(final int sample_index,final SAMSequenceDictionary dict,final WigItem item,int pos1) {
				this.sample_index = sample_index;
				this.contig = item.getChromosome();
				this.tid = dict.getSequenceIndex(this.contig);
				if(this.tid==-1) throw new JvarkitException.ContigNotFoundInDictionary(this.contig, dict);
				this.pos1 = pos1;
				this.rawdp = (int)item.getWigValue();
				}
			}
	
	private static class OneSampleDepthIterator
		extends AbstractIterator<OneSampleDepth>
		implements Closeable
		{
		final BufferedReader br ;
		final int sample_index;
		final SAMSequenceDictionary dict;
		OneSampleDepthIterator(final int sample_index,SAMSequenceDictionary dict,final File f) throws IOException {
			br= IOUtils.openFileForBufferedReading(f);
			this.sample_index = sample_index;
			this.dict = dict;
			}
		@Override
		protected OneSampleDepth advance()
			{
			try {
				final String line  = br.readLine();
				if(line==null) return null;
				return new OneSampleDepth(this.sample_index, this.dict, line);
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public void close() throws IOException
			{
			CloserUtil.close(br);
			}
		}
	
	/** iterate over the values of a bigwig/bigbed file */
	private static class OneSampleBigWigIterator
	extends AbstractIterator<OneSampleDepth>
	implements Closeable
		{
		private BBFileReader bbFileReader=null;
		final int sample_index;
		final SAMSequenceDictionary dict;
		final BigWigIterator bwIter;
		WigItem current_item = null;
		int position_in_current_item = -1; 
		OneSampleBigWigIterator(
				final int sample_index,
				SAMSequenceDictionary dict,
				final String path) throws IOException {
			this.bbFileReader = new BBFileReader(path);
			this.sample_index = sample_index;
			this.dict = dict;
			this.bwIter = this.bbFileReader.getBigWigIterator();
			}
		@Override
		protected OneSampleDepth advance()
			{
			if(current_item!=null)
				{
				if(position_in_current_item> current_item.getEndBase())
					{
					current_item = null;
					position_in_current_item = -1;
					}
				else
					{
					final OneSampleDepth osd = new OneSampleDepth(
							this.sample_index,
							this.dict,
							this.current_item,
							this.position_in_current_item
							);
					this.position_in_current_item++;
					return osd;
					}
				}
			try {
				if(this.bwIter.hasNext()) {
					close();
					return null;
					}
				this.current_item =this.bwIter.next();
				if(this.current_item==null) {
					close();
					return null;
					}
				this.position_in_current_item = this.current_item.getStartBase()+1;/* +1 for next iteration */
				return new OneSampleDepth(
						this.sample_index,
						this.dict,
						this.current_item,
						 this.current_item.getStartBase()
						);
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public void close() throws IOException
			{
			this.current_item=null;
			this.position_in_current_item=-1;
			CloserUtil.close(this.bwIter);
			CloserUtil.close(this.bbFileReader.getBBFis());
			CloserUtil.close(this.bbFileReader);
			}
		}

	
	private class MultipleSampleDepthIterators
		extends AbstractIterator<DepthLine>
		implements Closeable
		{
		final CharSplitter tab = CharSplitter.TAB;
		final BufferedReader br;
		MultipleSampleDepthIterators(final BufferedReader br) {
			this.br = br;
			}
		@Override
		protected DepthLine advance()
			{
			try {
				String line;
				for(;;) {
					line = br.readLine();
					if(line==null) return null;
					if(!line.startsWith("#")) break;
					}
				final String tokens[] = this.tab.split(line);
				if(tokens.length<3) {
					throw new JvarkitException.TokenErrors("expected at least 3 words",tokens);
					}
				final int n_samples = tokens.length-2;
				final String contig=tokens[0];
				final int pos1 = Integer.parseInt(tokens[1]);
				final DepthLine depthLine = new DepthLine(contig, pos1,n_samples);
				for(int x=2;x<tokens.length;++x)
					{
					final int rawdp =  Integer.parseInt(tokens[x]);
					depthLine.depths[x-2] = rawdp;
					}
				
				return depthLine;
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public void close() throws IOException
			{
			CloserUtil.close(br);
			}
		}
	
	private class CombineOneSampleDepthIterators
		extends AbstractIterator<DepthLine>
		implements Closeable
		{
		private List<Iterator<OneSampleDepth>> oneSampleDepthIterators;
		private MergingIterator<OneSampleDepth> mergingIter;
		private EqualRangeIterator<OneSampleDepth> equal_range;
		public CombineOneSampleDepthIterators(final SAMSequenceDictionary dic,final List<File> oneSampleFiles ) throws IOException
			{
			this.oneSampleDepthIterators = new ArrayList<>(oneSampleFiles.size());
			for(int i=0;i< oneSampleFiles.size();++i)
				{
				final File sampleFile=oneSampleFiles.get(i);
				final Iterator<OneSampleDepth> oiter;
				if(sampleFile.getName().toLowerCase().endsWith(".bw") ||
					sampleFile.getName().toLowerCase().endsWith(".bigwig") )
					{
					oiter = new OneSampleBigWigIterator(i, dic, sampleFile.getPath());
					}
				else
					{
					oiter = new OneSampleDepthIterator(i, dic,oneSampleFiles.get(i));
					}
				oneSampleDepthIterators.add(oiter);
				}
			final Comparator<OneSampleDepth> cmp= (A, B)->{
				final int d = Integer.compare(A.tid,B.tid);
				if(d!=0) return d;
				return Integer.compare(A.pos1,B.pos1);				
				}; 
			this.mergingIter = new MergingIterator<>(cmp,oneSampleDepthIterators);
			this.equal_range =  new EqualRangeIterator<>(this.mergingIter,cmp);
			}
		
		@Override
		protected DepthLine advance()
			{
			if(!this.equal_range.hasNext()) return null;
			final List<OneSampleDepth> row = this.equal_range.next();
			final DepthLine dp = new DepthLine(row.get(0).contig, row.get(0).pos1, oneSampleDepthIterators.size());
			
			for(final OneSampleDepth osd: row)
				{
				dp.depths[osd.sample_index] = osd.rawdp;
				}
			return dp;
			}
		@Override
		public void close() throws IOException
			{
			CloserUtil.close(this.equal_range);
			CloserUtil.close(this.mergingIter);
			CloserUtil.close(this.oneSampleDepthIterators);
			}
		}
	
	public static void main(final String[] args) {
		new NaiveCnvDetector().instanceMainWithExit(args);
		}
	}
