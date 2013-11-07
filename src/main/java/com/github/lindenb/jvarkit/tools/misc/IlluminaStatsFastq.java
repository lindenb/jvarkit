package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.regex.Pattern;

import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.samtools.SAMUtils;

import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.cli.GetOpt;
import com.github.lindenb.jvarkit.util.illumina.FastQName;

public class IlluminaStatsFastq
	extends AbstractCommandLineProgram
	{
	
	
	
	
	private static class Bases
		{
		long A=0L;
		long T=0L;
		long G=0L;
		long C=0L;
		long N=0L;
		}
	
	/* archive factory, where to put the results */
    private ArchiveFactory archiveFactory;
    private PrintWriter wnames=null;
    private PrintWriter wcount=null;
    private PrintWriter wquals=null;
    private PrintWriter wbadfastq=null;
    private PrintWriter whistquals=null;
    private PrintWriter wqualperpos=null;
    private PrintWriter wbases=null;
    private PrintWriter wlength=null;
    private PrintWriter wDNAIndexes=null;
    	
    private final Pattern DNARegex=Pattern.compile("[ATGCatcNn]{4,8}");
	
	private static void tsv(PrintWriter out,Object...array)
		{
		for(int i=0;i< array.length;++i)
			{
			if(i>0) out.print('\t');
			out.print(array[i]);
			}
		out.println();
		}
	
	@Override
	public String getProgramName()
		{
		return "Reads filenames from stdin: Count FASTQs in Illumina Result.";
		}
	
	private void analyze(File f) throws IOException
		{
		if(f==null) return;
		if(!(f.getName().endsWith(".fastq.gz") && f.isFile())) return;
		
			
			final int QUALITY_STEP=5;
			
			
			
			info(f.toString());
			FastQName fq=FastQName.parse(f);
			
			Counter<Integer> qualityHistogram=new Counter<Integer>();
			Counter<Integer> pos2quality=new Counter<Integer>();
			List<Bases> pos2bases=new ArrayList<Bases>(300);
			Counter<Integer> lengths=new Counter<Integer>();
			Counter<Integer> pos2count=new Counter<Integer>();
			Counter<String> dnaIndexes=new Counter<String>();
			long nReads=0L;
			double sum_qualities=0L;
			long count_bases=0L;
			long count_read_fails_filter=0L;
			long count_read_doesnt_fail_filter=0L;
			FastqReader r=null;
			try
				{
				r=new FastqReader(f);
				while(r.hasNext())
					{
					FastqRecord record=r.next();
					++nReads;
					if(record.getReadHeader().contains(":Y:"))
						{
						count_read_fails_filter++;
						continue;
						}
					else if(record.getReadHeader().contains(":N:"))
						{
						count_read_doesnt_fail_filter++;
						}
					
					if(COUNT_INDEX>0)
						{
						//index
						int last_colon=record.getReadHeader().lastIndexOf(':');
						if(last_colon!=-1 && last_colon+1< record.getReadHeader().length())
							{
							String dnaIndex=record.getReadHeader().substring(last_colon+1).trim().toUpperCase();
							if(DNARegex.matcher(dnaIndex).matches())
								{
								dnaIndexes.incr(dnaIndex);
								}
							}
						}
					
					byte phred[]=SAMUtils.fastqToPhred(record.getBaseQualityString());
					
					for(int i=0;i< phred.length ;++i)
						{
						sum_qualities+=phred[i];
						count_bases++;
						
						qualityHistogram.incr(phred[i]/QUALITY_STEP);
						pos2quality.incr(i,phred[i]);
						pos2count.incr(i);
						}
					/* get base usage */
					while(pos2bases.size() <record.getReadString().length())
						{
						pos2bases.add(new Bases());
						}
					for(int i=0;i< record.getReadString().length() ;++i)
						{
						Bases bases=pos2bases.get(i);
						switch(record.getReadString().charAt(i))
							{
							case 'A': case 'a':bases.A++;break;
							case 'T': case 't':bases.T++;break;
							case 'G': case 'g':bases.G++;break;
							case 'C': case 'c':bases.C++;break;
							default: bases.N++;break;
							}
						}
					lengths.incr(record.getBaseQualityString().length());
					}
				}
			catch(Exception err2)
				{
				error(err2,"BOUM "+err2.getMessage());
				err2.printStackTrace();
				tsv(wbadfastq,f.getPath(),err2.getMessage());
				return;
				}	
			finally
				{
				if(r!=null) r.close();
				r=null;
				}
			
			
			if(fq.isValid())
				{
				tsv(this.wnames,
					f.getPath(),
					(fq.isUndetermined()?"Undetermined":fq.getSample()),
					fq.getSeqIndex(),
					fq.getLane(),
					fq.getSide(),
					fq.getSplit(),
					fq.getFile().length()
					);
				}
			else
				{
				tsv(wbadfastq,f.getPath());
				}
			
			tsv(this.wcount,
				f.getPath(),
				nReads,
				count_read_fails_filter,
				count_read_doesnt_fail_filter
				);
			
			tsv(this.wquals,
				f.getPath(),
				sum_qualities/count_bases
				);
			for(Integer step:qualityHistogram.keySet())
				{
				tsv(this.whistquals,
						f.getPath(),
						step*QUALITY_STEP,
						qualityHistogram.count(step)
						);
				
				}
			for(Integer position:pos2quality.keySet())
				{
				tsv(this.wqualperpos,
						f.getPath(),
						position+1,
						pos2quality.count(position)/(double)pos2count.count(position),
						pos2count.count(position)
						);
				}
			
			for(int i=0;i< pos2bases.size();++i)
				{
				Bases b=pos2bases.get(i);
				tsv(this.wbases,
					f.getPath(),
					i+1,b.A,b.T,b.G,b.C,b.N
					);
				}
			for(Integer L:lengths.keySet())
				{
				tsv(this.wlength,
						f.getPath(),
						L,
						lengths.count(L)
						);
				}
			
			int count_out=0;
			for(String dna:dnaIndexes.keySetDecreasing())
				{
				if(++count_out>COUNT_INDEX) break;
				tsv(this.wDNAIndexes,f.getPath(),dna,dnaIndexes.count(dna));
				}
			
			
		}

	private int COUNT_INDEX=0;
	private File OUT=null;
	
	@Override
	public void printUsage(PrintStream out)
		{
		printStandardPreamble(out);
		out.println("Options:");
		out.println(" -h help. This Screen");
		out.println(" -v print version and exits");
		out.println(" -L (log-level , a "+Level.class.getName()+"). Optional");
		out.println(" -X (int) maximum number of DNA indexes to print. memory consuming if not 0. Optional");
		out.println(" -o (filename out) Directory or ZIP. Required.");
		}
	
	@Override
	public int doWork(String args[])
		{
		GetOpt getopt=new GetOpt();
		int c;
		while((c=getopt.getopt(args, "hvL:o:X:"))!=-1)
			{
			switch(c)
				{
				case 'h': printUsage();return 0;
				case 'v': System.out.println(getVersion());break;
				case 'L': getLogger().setLevel(Level.parse(getopt.getOptArg()));break;
				case 'X':
					{
					COUNT_INDEX=Integer.parseInt(getopt.getOptArg());
					break;
					}
				case 'o':
					{
					this.OUT=new File(getopt.getOptArg());
					break;
					}
				default:
					{
					System.err.println("Unknown option or missing argument. "+(char)getopt.getOptOpt());
					return -1;
					}
				}
			}
		
		
		if(getopt.getOptInd()!=args.length)
			{
			System.err.println("Expected reads from stdin. Illegal Number of arguments.");
			return -1;
			}
		if(OUT==null)
			{
			System.err.println("undefined output file.");
			return -1;
			}
		
		
		try {
			
			archiveFactory=ArchiveFactory.open(OUT);
			this.wnames = archiveFactory.openWriter("names.tsv");
			this.wcount = archiveFactory.openWriter("counts.tsv");
			this.wquals = archiveFactory.openWriter("quals.tsv");
			this.wbadfastq = archiveFactory.openWriter("notfastq.tsv");
			this.whistquals = archiveFactory.openWriter("histquals.tsv");
			this.wqualperpos = archiveFactory.openWriter("histpos2qual.tsv");
			this.wbases = archiveFactory.openWriter("bases.tsv");
			this.wlength = archiveFactory.openWriter("lengths.tsv");
			this.wDNAIndexes = archiveFactory.openWriter("indexes.tsv");
			
			info("reading from stdin");
			BufferedReader in=new BufferedReader(new InputStreamReader(System.in));
			String line;
			while((line=in.readLine())!=null)
				{	
				if(line.isEmpty() || line.startsWith("#")) continue;
				try 
					{
					analyze(new File(line));
					}
				catch(IOException err)
					{
					warning(err,"Cannot analyse "+line);
					}
				}
			
			in.close();
			
			
			for(PrintWriter pw: new PrintWriter[]{
					this.wnames,
					this.wcount,
					this.wquals,
					this.wbadfastq,
					this.whistquals,
					this.wqualperpos,
					this.wbases,
					this.wlength,
					this.wDNAIndexes})
				{
				pw.flush();
				pw.close();
				}
			

			
			if(archiveFactory!=null)
				{
				this.archiveFactory.close();
				}
			} 
		catch (Exception e) {
			e.printStackTrace();
			error(e,"ERROR:"+e.getMessage());
			return -1;
			}
		finally	
			{
			
			}
		return 0;
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new IlluminaStatsFastq().instanceMainWithExit(args);
		}

}
