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
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.inference.ChiSquareTest;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.tools.vcfviewgui.PedFile.Sample;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;


/**
BEGIN_DOC

Input is the output of samtools depth.

END_DOC

 */
@Program(name="naivecnvdetector",
	description="experimental CNV detection for multiple samples. Doesn't work for now.",
	keywords= {"cnv","bam","sam"},
	generate_doc=false
	)
public class NaiveCnvDetector extends Launcher
	{
	private static final Logger LOG = Logger.build(NaiveCnvDetector.class).make();
	
	//@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	//private File outputFile=null;
	
	@Parameter(names={"-c"},description="config file. Tab delimited. Sample-name(tab)mean-depth(tab)integer[affected=1,non-affected=0]",required=true)
	private File configFile=null;	
	/** size of a window */
	@Parameter(names={"-w"},description="window size")
	private int windowSize=1000;
	@Parameter(names={"-s"},description="window shift")
	private int windowShift=500;

	@Parameter(names={"--weirdDepth"},description="Treat normalized depth greater than this value as 'weird' and discard the sliding windows at this place.")
	private int weirdDepth=500;
	@Parameter(names={"-nafd","--no-affected-for-depth"},description="When calculating the average depth in one window, do not include the affected samples. use when the number of affected ~= non-affected.")
	private boolean non_affected_for_depth=false;
	@Parameter(names={"-md","--min-dp"},description="At least one sample must have a normalized-depth greater than this value.")
	private int min_depth = 20;
	
	private final class SampleInfo
		{
		String name;
		int index;
		double meanDepth=0.0;
		double adjustDepth = 1.0;
		boolean affected=false;
		}
	
	private final List<SampleInfo> sampleList  = new ArrayList<>();
	
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
	
	private void dump() {
		if(depthBuffer.isEmpty()) return;
		DepthInterval rec = new DepthInterval(this.depthBuffer);
		if(Arrays.stream(rec.depths).noneMatch(V->V>=min_depth)) return;
		if(Arrays.stream(rec.depths).anyMatch(V->V>=weirdDepth)) return;
		if(Arrays.stream(rec.depths,0,6).max().getAsDouble()<25) return;
		
		
		final StandardDeviation standardDeviation=new StandardDeviation();
		Arrays.stream(rec.depths).forEach(V->standardDeviation.increment(V));
		standardDeviation.evaluate();
		
		final double median_depth = new Median().evaluate(rec.depths,0,6);
		final double dup = median_depth * 2.8;
		final double del = median_depth * 0.2;
		List<Integer> dupIndexes = new ArrayList<>();
		List<Integer> delIndexes = new ArrayList<>();
		for(int i=6;i< rec.depths.length;i++)
			{
			double d = rec.depths[i];
			if(d>dup)
				{
				dupIndexes.add(i);
				}
			if(d<del)
				{
				delIndexes.add(i);
				}
			}
		String msg ;
		if(delIndexes.size()>=1 && dupIndexes.isEmpty())
			{
			msg=(rec.contig+"\t"+rec.start+"\t"+rec.end+"\tDEL\t["+delIndexes.get(0)+"]");
			}
		else if(dupIndexes.size()>=1 && delIndexes.isEmpty())
			{
			msg=(rec.contig+"\t"+rec.start+"\t"+rec.end+"\tDUP\t["+dupIndexes.get(0)+"]");
			}
		else
			{
			return;
			}
		msg+="\t"+Arrays.stream(rec.depths).boxed().map(V->V.toString()).collect(Collectors.joining(" , "));
		
		final ChiSquareTest chiSquareTest = new ChiSquareTest();
		final double p_value=chiSquareTest.chiSquare(new long[][]{{0L},{}});
		System.out.println(msg);
		
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
		int num_samples = -1;
		BufferedReader samDepthReader=null;
		
		try
			{
			final CharSplitter tab = CharSplitter.TAB;
			
			
			this.sampleList.addAll( IOUtil.slurpLines(this.configFile).stream().
					filter(S->!(StringUtil.isBlank(S) || S.startsWith("#"))).
					map(S->CharSplitter.TAB.split(S)).
					map(T->{
						final SampleInfo si=new SampleInfo();
						si.name = T[0];
						si.meanDepth = Double.parseDouble(T[1]);
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
				this.sampleList.get(i).adjustDepth=this.sampleList.get(i).meanDepth* max_depth;
				}
			
			
			samDepthReader = super.openBufferedReader(oneFileOrNull(args));
			String line;
			
			
			
			while((line=samDepthReader.readLine())!=null)
				{
				final String tokens[] = tab.split(line);
				if(tokens.length<3) {
					throw new JvarkitException.TokenErrors("expected at least 3 words",tokens);
					}
				if(this.sampleList.size()+2!=tokens.length)
					{
					throw new JvarkitException.TokenErrors("expected at least "+(num_samples+3)+" words",tokens);
					}
				final String contig=tokens[0];
				final int pos1 = Integer.parseInt(tokens[1]);
				final DepthLine depthLine = new DepthLine(contig, pos1, num_samples);
				for(int x=2;x<tokens.length;++x)
					{
					depthLine.depths[x-2] = this.sampleList.get(x-2).adjustDepth * Integer.parseInt(tokens[x]);
					}
				if(!this.depthBuffer.isEmpty())
					{
					final DepthLine last = this.depthBuffer.get(this.depthBuffer.size()-1);
					if(!last.contig.equals(depthLine.contig))
						{
						dump();
						this.depthBuffer.clear();
						}
					else if(last.pos+1!=depthLine.pos)
						{
						dump();
						this.depthBuffer.clear();
						}
					else if(last.pos>=depthLine.pos)
						{
						dump();
						this.depthBuffer.clear();
						}
					}
				this.depthBuffer.add(depthLine);
				if(this.depthBuffer.size()==this.windowSize)
					{
					dump();
					for(int x=0;x<this.windowShift && !this.depthBuffer.isEmpty();++x)
						{
						this.depthBuffer.remove(0);
						}
					}
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
			CloserUtil.close(samDepthReader);
			}
		}
	
	public static void main(final String[] args) {
		new NaiveCnvDetector().instanceMainWithExit(args);
		}
	}
