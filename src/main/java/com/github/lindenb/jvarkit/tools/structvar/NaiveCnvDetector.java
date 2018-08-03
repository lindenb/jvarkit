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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.BufferedReader;
import java.io.File;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.inference.ChiSquareTest;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;


/**
BEGIN_DOC

Input is the output of samtools depth.


## Example

```
samtools depth -r '1:1234-567' *.bam |\
	java -jar dist/naivecnvdetector.jar  > out.tsv
```


END_DOC

 */
@Program(name="naivecnvdetector",
	description="experimental CNV detection for multiple samples.",
	keywords= {"cnv","bam","sam"}
	)
public class NaiveCnvDetector extends Launcher
	{
	private static final Logger LOG = Logger.build(NaiveCnvDetector.class).make();
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@Parameter(names={"-c"},description="config file. Tab delimited. Sample-name(tab)mean-depth(tab)integer[affected=1,non-affected=0]. If this file is not specified , all samples are considered unaffected (discovery mode).")
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
	@Parameter(names={"-del","--del","--deletion"},description="Deletion Treshold. Which fraction of the median depth is considered as aa deletion. Must be <1.0" )
	private double deletion_treshold = 0.5;
	@Parameter(names={"-dup","--dup","--duplication"},description="Duplication Treshold. Which fraction of the median depth is considered as a duplication. Must be >1.0" )
	private double duplication_treshold = 1.9;
	@Parameter(names={"-disable-both"},description="Disable the following criteria: there cannot be a DEL and a DUP at the same place." )
	private boolean disable_both_del_dup = false;


	
	private class SampleInfo
		{
		String name;
		int index;
		double meanDepth=0.0;
		double adjustDepth = 1.0;
		boolean affected=false;
		
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
		
		if(median_unaffected_depth<0) return;
		
		final List<SampleInfo> delSamples = this.sampleList.
				stream().
				filter(SI-> rec.depths[SI.index] < median_unaffected_depth*this.deletion_treshold).
				collect(Collectors.toList())
				;
		final List<SampleInfo> dupSamples = this.sampleList.
				stream().
				filter(SI-> rec.depths[SI.index] > median_unaffected_depth*this.duplication_treshold).
				collect(Collectors.toList())
				;
		
		if(delSamples.isEmpty() && dupSamples.isEmpty()) return;
		
		final Set<SampleInfo> noCnvSamples = new HashSet<>(this.sampleList);
		noCnvSamples.removeAll(delSamples);
		noCnvSamples.removeAll(dupSamples);
		
		
		//interval contains DEL *and* DUP
		if(!this.disable_both_del_dup && 
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
		out.print("\t");
		
		
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
			out.print(p_value);
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
			String prefix=side==0?"DEL":"DUP";
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
		out.print("\t");
		
		if(this.count_affected_samples>0) {
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
		
		BufferedReader samDepthReader=null;
		PrintWriter out = null;
		try
			{
			
			final CharSplitter tab = CharSplitter.TAB;
			
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
			
			
			samDepthReader = super.openBufferedReader(oneFileOrNull(args));
			String line;
			
			
			
			while((line=samDepthReader.readLine())!=null)
				{
				final String tokens[] = tab.split(line);
				if(tokens.length<3) {
					throw new JvarkitException.TokenErrors("expected at least 3 words",tokens);
					}
				if(this.sampleList.isEmpty())
					{
					LOG.info("building default sample list");
					this.count_affected_samples = 0;
					this.count_unaffected_samples = tokens.length-2;
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
				else if(this.sampleList.size()+2!=tokens.length)
					{
					throw new JvarkitException.TokenErrors("expected at least "+( this.sampleList.size()+3)+" words",tokens);
					}
				final String contig=tokens[0];
				final int pos1 = Integer.parseInt(tokens[1]);
				final DepthLine depthLine = new DepthLine(contig, pos1, this.sampleList.size());
				for(int x=2;x<tokens.length;++x)
					{
					depthLine.depths[x-2] = this.sampleList.get(x-2).adjustDepth * Integer.parseInt(tokens[x]);
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
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(samDepthReader);
			CloserUtil.close(out);
			}
		}
	
	public static void main(final String[] args) {
		new NaiveCnvDetector().instanceMainWithExit(args);
		}
	}
