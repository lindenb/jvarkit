package com.github.lindenb.jvarkit.tools.biostar;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
BEGIN_DOC

## input

input is a bed file

## see also

  * https://gist.github.com/lindenb/6130880/

END_DOC
 */
@Program(name="biostar77828",
description="Divide the human genome among X cores, taking into account gaps",
		biostars= {77828,369434},
		keywords={"workflow","reference","parallel"}
		)
public class Biostar77828 extends Launcher
	{

	private static final Logger LOG = Logger.build(Biostar77828.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-minc","--minc"},description="min core")
	private int MINC = 20 ;
	@Parameter(names={"-maxc","--maxc"},description="max core")
	private int MAXC = 30 ;
	@Parameter(names={"-iter","--iter"},description="number of iterations")
	private long N_ITERATIONS = 1000000 ;



    
    private List<Segment> all_segments=new ArrayList<Segment>();
    private final Random random=new Random(System.currentTimeMillis());
    private long effective_genome_size=0L;
    	
    private static class Segment implements Comparable<Segment>
    	{
    	final String chrom;
    	final int start;
    	final int end;
    	final String name;
    	
    	Segment(String chrom,int start,int end)
    		{
    		this(chrom,start,end,chrom+":"+start+":"+end);
    		}
    	
    	Segment(String chrom,int start,int end,String name)
			{
			this.chrom=chrom;
			this.start=start;
			this.end=end;
			this.name=name;
			}
    	
    	int size()
    		{
    		return end-start;
    		}
    	@Override
    	public int compareTo(final Segment o)
    		{
    		int i=chrom.compareTo(o.chrom);
    		if(i!=0) return i;
    		i=start-o.start;
    		if(i!=0) return i;
    		return end-o.end;
    		}
    	}
    
    private static class Core
    	{
    	final List<Segment> segments=new ArrayList<Segment>();
    	long length()
    		{
    		return this.segments.stream().mapToLong(S->S.size()).sum();
    		}
    	void simplify()
			{
			Collections.sort(this.segments);
			int i=0;
			while(i+1< this.segments.size())
				{
				final Segment S1=this.segments.get(i);
				final Segment S2=this.segments.get(i+1);
				if(S1.chrom.equals(S2.chrom) && S1.end==S2.start)
					{
					this.segments.set(i, new Segment(S1.chrom, S1.start, S1.end,S1.name));
					this.segments.remove(i+1);
					}
				else
					{
					i++;
					}
				}
			}
    	}
    
    private static class Solution
    	implements Comparable<Solution>
    	{
    	long generation=0L;
    	List<Core> cores=new ArrayList<Core>();
    	private Double _mean=null;
    	private Double _stddev=null;
    	double mean()
			{
			if(_mean==null)
				{
				_mean=0.0;
				for(Core s:this.cores) _mean+=s.length();
				_mean/=this.cores.size();
				}
			return _mean;
			}

    	Double stddev()
    		{
    		if(_stddev==null)
				{
    			double m=mean();
    			_stddev=0.0;
				for(Core s:this.cores) _stddev+=Math.pow(s.length()-m,2.0);
				_stddev=Math.sqrt(_stddev/(this.cores.size()));
				}
			return _stddev;    		
			}
    	
    	@Override
    	public int compareTo(Solution other)
    		{
    		return this.stddev().compareTo(other.stddev());
    		}
    	
    
    	
    	public void print(PrintStream out)
    		{
    		String header="##mean:"+mean()+"\n##stdev:"+stddev()+" \n##cores:"+this.cores.size()+"\n##generation:"+this.generation;
    		out.println(header);
    		for(int i=0;i< this.cores.size();++i)
    			{
    			this.cores.get(i).simplify();
    			
    			out.println("#cores["+i+"]. N="+this.cores.get(i).segments.size()+"  size_bp="+this.cores.get(i).length()+" bp.");
    			for(Segment seg:this.cores.get(i).segments)
    				{
    				out.print(seg.chrom);
    				out.print('\t');
    				out.print(seg.start);
    				out.print('\t');
    				out.print(seg.end);
    				out.print('\t');
    				out.print(seg.size());
    				out.print('\t');
    				out.print(seg.name);
    				out.println();
    				}
    			}
    		out.println(header);
    		}
    	@Override
    	public String toString() {
    		return "cores:"+this.cores.size()+" mean "+mean()+" stddev:"+stddev()+" gen:"+generation;
    		}
    	}
    
    private Solution createSolution()
    	{
    	int n_cores=
    			this.MINC+(this.MINC>=this.MAXC?0:this.random.nextInt(this.MAXC-this.MINC))
    			;
    	long optimal_size= (this.effective_genome_size/n_cores);
    	
    	Solution sol=new Solution();
    	
    	//create a copy, split the segments
    	ArrayList<Segment> segments=new ArrayList<Segment>(this.all_segments.size());
    	for(Segment seg:this.all_segments)
    		{
    		if((long)seg.size()<=optimal_size)
    			{
    			segments.add(seg);
    			}
    		else
    			{
    			int split_count=(int)(Math.ceil(seg.size()/(double)optimal_size));
    			int fragment_size=seg.size()/split_count;
    			int start=seg.start;
    			for(int s=0;s<split_count;++s)
    				{
    				int end=(s+1==split_count?seg.end:start+fragment_size);
    				segments.add(new Segment(seg.chrom, start, end,seg.name));
    				start=end;
    				}
    			}
    		}
    	//shuffle
    	Collections.shuffle(segments, this.random);
    	
    	//create cores, and add one segment 
    	for(int i=0;i< n_cores && i< segments.size();++i)
    		{
    		//get last
    		Core core=new Core();
    		core.segments.add(segments.get(i));
    		sol.cores.add(core);
    		}
    	
    	double seq_total_size=0.0;
    	int seq_count=0;
    	//fill core with random segments
    	for(int i=sol.cores.size();i< segments.size();++i)
    		{
    		Segment seg=segments.get(i);
    		//start a new sequence of insertion
    		if(seq_count%sol.cores.size()==0)
    			{
        		//put it in a random core
        		Core core=sol.cores.get(this.random.nextInt(sol.cores.size()));
        		core.segments.add(seg);
        		seq_total_size=seg.size();
        		seq_count=1;
    			}
    		else
    			{
    			double mean=seq_total_size/seq_count;
    			//find the best core to put this segment best=smaller distance to mean
    			Core best=null;
    			for(Core core:sol.cores)
    				{
    				if(best==null ||
    					(Math.abs((core.length()+seg.size())-mean) < Math.abs((best.length()+seg.size())-mean)))
    					{
    					best=core;
    					}
    			
    				}
    			best.segments.add(seg);
    			seq_total_size+=seg.size();
        		seq_count++;
    			}
    		
    		}
    	
    	return sol;
    	}
    @Override
    public int doWork(final List<String> args) {
		PrintStream pw =null;
    	try
	    	{
    		
	    	LOG.info("load BED");
	    	BufferedReader in=super.openBufferedReader(oneFileOrNull(args));
	    
	    	final BedLineCodec codec = new BedLineCodec();
	    	String line;
	    	while((line=in.readLine())!=null)
				{		
	    		if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
	    		final BedLine bedLine = codec.decode(line);
	    		if(bedLine==null) continue;
    			if(bedLine.getColumnCount()<3) throw new IOException("bad BED input "+bedLine);
    			final Segment seg=new Segment(
    					bedLine.getContig(),
    					bedLine.getStart()-1,
    					bedLine.getEnd()
    					)
    					;
    			if(seg.size()==0) continue;
    			this.effective_genome_size+=seg.size();
    			this.all_segments.add(seg);
				}
	    	
	    	pw = super.openFileOrStdoutAsPrintStream(this.outputFile);
	    	
	    	Solution best=null;
	    	for(long generation=0;generation< this.N_ITERATIONS;++generation)
	    		{
	    		if(generation%100000==0) LOG.info("generation:"+generation+"/"+this.N_ITERATIONS+" "+
	    				(int)((generation/(double)this.N_ITERATIONS)*100.0)+"% "+String.valueOf(best));
	    		Solution sol=createSolution();
	    		sol.generation=generation;
	    		if(best==null || sol.compareTo(best)<0)
	    			{
	    			best=sol;
	    			/*
	    			if(LOG.isDebugEnabled())
	    				{
	    				LOG.info("%%generation:"+generation);
	    				best.print(stderr());
	    				}*/
	    			}
	    		}
	    	if(best!=null)
	    		{
	    		best.print(pw);
	    		}
	    	pw.flush();
	    	pw.close();
	    	return RETURN_OK;
	    	}
    	catch(final Exception err)
    		{
    		LOG.error(err);
    		return -1;
    		}
    	finally {
    		CloserUtil.close(pw);
    	}
    }

    
	public static void main(final String[] args) {
		new Biostar77828().instanceMainWithExit(args);

	}

}
