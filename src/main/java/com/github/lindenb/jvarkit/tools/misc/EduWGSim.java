package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import net.sf.picard.fastq.BasicFastqWriter;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;
import net.sf.picard.util.Interval;

import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;



public class EduWGSim extends AbstractCommandLineProgram
	{
	public int READ_LENGTH=100; 
	public int FRAGMENT_SIZE=500;
	public int N=10000; 
	
	private static class Mut
		{
		String chrom;
		int pos1;
		char base;
		}
	private Random random=new Random();
	
	
	@Override
	protected int doWork()
		{
		File f1=null;
		FastqWriter fq1=new BasicFastqWriter(System.out);
		Mut m=new Mut();
		List<Interval> intervals=new ArrayList<Interval>();
		int n=0;
		while(n<N)
			{
			Interval interval=intervals.get(this.random.nextInt(intervals.size()));
			int pos=this.random.nextInt(interval.length()-FRAGMENT_SIZE);
			
			++n;
			}
		fq1.close();
		return 0;
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new EduWGSim().instanceMainWithExit(args);
		}

}
