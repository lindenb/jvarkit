package com.github.lindenb.jvarkit.tools.biostar;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.regex.Pattern;

import htsjdk.samtools.cmdline.Option;
import htsjdk.samtools.cmdline.StandardOptionDefinitions;
import htsjdk.samtools.cmdline.Usage;
import htsjdk.samtools.util.Log;

import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;

public class Biostar77828 extends AbstractCommandLineProgram
	{
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+
		" Divide the human genome among X cores, taking into account gaps See http://www.biostars.org/p/77828/ .";
	private static final Log LOG=Log.getInstance(Biostar77828.class);



    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME,doc="Bed file input (default stdin)",optional=true)
    public File IN=null;
    
    @Option(shortName="MINC",doc="min_core",optional=true)
    public int MIN_CORE=20;
    @Option(shortName="MAXC",doc="max_core",optional=true)
    public int MAX_CORE=30;

    @Option(shortName="ITER",doc="number of iterations",optional=true)
    public long N_ITERATIONS=1000000;

    
    
    private List<Segment> all_segments=new ArrayList<Segment>();
    private final Random random=new Random(System.currentTimeMillis());
    private long effective_genome_size=0L;
    	
    private static class Segment implements Comparable<Segment>
    	{
    	String chrom;
    	int start;
    	int end;
    	String name;
    	
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
    	public int compareTo(Segment o)
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
    	List<Segment> segments=new ArrayList<Segment>();
    	long length()
    		{
    		long L=0L;
    		for(Segment s:this.segments) L+=s.size();
    		return L;
    		}
    	void simplify()
			{
			Collections.sort(this.segments);
			int i=0;
			while(i+1< this.segments.size())
				{
				Segment S1=this.segments.get(i);
				Segment S2=this.segments.get(i+1);
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
    			MIN_CORE+(MIN_CORE>=MAX_CORE?0:this.random.nextInt(MAX_CORE-MIN_CORE))
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
    protected int doWork()
    	{
    	try
	    	{
	    	LOG.info("load BED");
	    	BufferedReader in=null;
	    	if(IN==null)
	    		{
	    		LOG.info("read stdin");
	    		in=new BufferedReader(new InputStreamReader(System.in));
	    		}
	    	else
	    		{
	    		LOG.info("read "+IN);
	    		in=new BufferedReader(new FileReader(IN));
	    		}
	    	Pattern tab=Pattern.compile("[\t]");
	    	String line;
	    	while((line=in.readLine())!=null)
				{		
	    		if(line.isEmpty() || line.startsWith("#")) continue;
	    		String tokens[]=tab.split(line,5);
    			if(tokens.length<3) throw new IOException("bad BED input "+Arrays.asList(tokens));
    			Segment seg=new Segment(
    					tokens[0],
    					Integer.parseInt(tokens[1]),
    					Integer.parseInt(tokens[2]))
    					;
    			if(seg.size()==0) continue;
    			this.effective_genome_size+=seg.size();
    			this.all_segments.add(seg);
				}
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
	    			if(super.VERBOSITY==Log.LogLevel.DEBUG)
	    				{
	    				System.err.println("%%generation:"+generation);
	    				best.print(System.err);
	    				}
	    			}
	    		}
	    	if(best!=null)
	    		{
	    		best.print(System.out);
	    		}
	    	
	    	}
    	catch(Exception err)
    		{
    		LOG.error(err);
    		return -1;
    		}
    	return 0;
    	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Biostar77828().instanceMainWithExit(args);

	}

}
