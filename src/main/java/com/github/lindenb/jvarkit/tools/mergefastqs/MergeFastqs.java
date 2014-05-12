package com.github.lindenb.jvarkit.tools.mergefastqs;

import java.awt.Point;
import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;

import com.github.lindenb.jvarkit.util.picard.cmdline.CommandLineProgram;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.Log;

public class MergeFastqs extends CommandLineProgram
	{
	private Log LOG=Log.getInstance(MergeFastqs.class);
	public File IN1=null;
	public File IN2=null;
	public int MATCH=10;
	
	
	
	private char complement(char c)
		{
        switch(c)
	        {
	        case 'a':case 'A':return 'T';
	        case 't':case 'T':return 'A';
	        case 'g':case 'G':return 'C';
	        case 'c':case 'C':return 'G';
	        default: return 'N';
	        }
		}
	
	@Override
	protected int doWork()
		{
		FastqReader fqr1=null;
		FastqReader fqr2=null;
		FastqWriter fw1=null;
		FastqWriter fw2=null;
		FastqWriter fw3=null;
		int matrix[]=new int[0];
		try
			{
			Point point=new Point();
			fqr1=new FastqReader(IN1);
			fqr2=new FastqReader(IN2);
			while(fqr1.hasNext())
				{
				if(!fqr2.hasNext()) throw new IllegalStateException();
				FastqRecord fq1=fqr1.next();
				FastqRecord fq2=fqr1.next();
				String X=fq1.getReadString();
				String Y=fq2.getReadString();
				
				int width=X.length()+1;
				int height=Y.length()+1;
				if(matrix.length< (width*height))
					{
					matrix=new int[width*height];
					}
				Arrays.fill(matrix, 0);
				int best=0;
				
				for(int x=1;x<=X.length();++x)
					{
                    char c1= X.charAt((x-1));
					for(int y=1;y<=Y.length();++y)
						{
						char c2= complement(Y.charAt(Y.length()-1-(y-1)));
						int idx=matrix[y*width+x];
						if(c1==c2)
							{
							matrix[idx]=matrix[(y-1)*width+(x-1)]+1;
							if(matrix[idx]>best)
								{
								best=matrix[idx];
								point.x=x;
								point.y=y;
								}
							}
						else
							{
							matrix[idx]=0;
							}
						}
					
					}


	            if(best< this.MATCH)
	            	{
	            	if(fw1!=null) fw1.write(fq1);
	            	if(fw2!=null) fw2.write(fq2);
	            	}
	            else
	            	{
	            	
	            	
	            	if(fw3!=null)
	            		{
	            		StringBuilder dna=new StringBuilder();
	            		StringBuilder qual=new StringBuilder();
	            		
	            		for(int x=0;x<point.x;++x)
	            			{
	            			dna.append(X.charAt(x));
	            			}
	            		for(int x=0;x< best;++x)
	            			{
	            			if((int)fq1.getBaseQualityString().charAt(0)< (int)fq2.getBaseQualityString().charAt(0))
	            				{
	            				dna.append(X.charAt(x));
	            				qual.append("TODO");
	            				}
	            			else
	            				{
	            				dna.append(X.charAt(x));
	            				qual.append("TODO");
	            				}
	            			}
	            		
	            		for(int y=point.y;y>=0;--y)
	            			{
	            			dna.append(Y.charAt(y));
	            			}
	            		
	            		FastqRecord merged=new FastqRecord(
	            				fq1.getReadHeader(),
	            				dna.toString(),
	            				fq1.getBaseQualityHeader(),
	            				qual.toString()
	            				);
	            		fw3.write(merged);
	            		}
	            	}				
				}
			if(fqr2.hasNext()) throw new IllegalStateException();
			return 0;
			}
		catch (Exception e)
			{
			return -1;
			}
		finally
			{
			if(fqr1!=null) fqr1.close();
			if(fqr2!=null) fqr2.close();
			if(fw1!=null) fw1.close();
			if(fw2!=null) fw2.close();
			if(fw3!=null) fw3.close();
			}
		}

	}
