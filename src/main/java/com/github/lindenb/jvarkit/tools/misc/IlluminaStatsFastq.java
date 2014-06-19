package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.regex.Pattern;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMUtils;

import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.cli.GetOpt;
import com.github.lindenb.jvarkit.util.illumina.FastQName;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;

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
    private PrintWriter wsqlite=null;
    private int nThreads=1;
    	
	
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
	public String getProgramDescription() {
		return "Reads filenames from stdin: Count FASTQs in Illumina Result.";
		}
	
	private class Analyzer extends Thread
		{
		private File fastqFile;
	    private final Pattern DNARegex=Pattern.compile("[ATGCatcNn]{4,8}");
	    private String hash;
	    
		Analyzer(File fastqFile)
			{
			this.fastqFile=fastqFile;
			 try {
				 MessageDigest md5=MessageDigest.getInstance("MD5");
		   		md5.reset();
		   		md5.update(String.valueOf(fastqFile.getPath()).getBytes());
		   	
		        
		   		hash = new BigInteger(1, md5.digest()).toString(16);
		        if (hash.length() != 32) {
		            final String zeros = "00000000000000000000000000000000";
		            hash= zeros.substring(0, 32 - hash.length()) + hash;
		        }
		         } catch (NoSuchAlgorithmException e) {
		           throw new RuntimeException("MD5 algorithm not found",e);
		         }
			}
		
		
			
		
   
		
		@Override
		public void run()
			{
			try {
				analyze(this.fastqFile);
				}
			catch (Exception e) {
				error(e);
				}
		
			}
		
		IlluminaStatsFastq owner()
			{
			return IlluminaStatsFastq.this;
			}
		
		private void analyze(File f) throws IOException
			{
			if(f==null) return;
			if(!(f.getName().endsWith(".fastq.gz") && f.isFile())) return;
			if(!f.canRead())
				{
				synchronized (IlluminaStatsFastq.class)
					{
					tsv(owner().wbadfastq,f.getPath(),"Cannot read");
					}
				return;
				}
				
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
					r=new FourLinesFastqReader(f);
					r.setValidationStringency(ValidationStringency.LENIENT);
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
						
						if(owner().COUNT_INDEX>0)
							{
							//index
							int last_colon=record.getReadHeader().lastIndexOf(':');
							if(last_colon!=-1 && last_colon+1< record.getReadHeader().length())
								{
								String dnaIndex=record.getReadHeader().substring(last_colon+1).trim().toUpperCase();
								if(this.DNARegex.matcher(dnaIndex).matches())
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
					synchronized (IlluminaStatsFastq.class)
						{
						tsv(owner().wbadfastq,f.getPath(),this.hash,err2.getMessage());
						}
					
					return;
					}	
				finally
					{
					if(r!=null) r.close();
					r=null;
					}
				
				synchronized (IlluminaStatsFastq.class)
					{
					if(fq.isValid())
						{
						tsv(owner().wnames,
							f.getPath(),
							f.getParentFile(),
							f.getName(),
							this.hash,
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
						tsv(owner().wbadfastq,f.getPath(),this.hash);
						}
					
					tsv(owner().wcount,
						this.hash,
						nReads,
						count_read_fails_filter,
						count_read_doesnt_fail_filter
						);
					
					tsv(owner().wquals,
						this.hash,
						sum_qualities/count_bases
						);
					for(Integer step:qualityHistogram.keySet())
						{
						tsv(owner().whistquals,
								this.hash,
								step*QUALITY_STEP,
								qualityHistogram.count(step)
								);
						
						}
					for(Integer position:pos2quality.keySet())
						{
						tsv(owner().wqualperpos,
								this.hash,
								position+1,
								pos2quality.count(position)/(double)pos2count.count(position),
								pos2count.count(position)
								);
						}
					
					for(int i=0;i< pos2bases.size();++i)
						{
						Bases b=pos2bases.get(i);
						tsv(owner().wbases,
							this.hash,
							i+1,b.A,b.T,b.G,b.C,b.N
							);
						}
					for(Integer L:lengths.keySet())
						{
						tsv(owner().wlength,
								this.hash,
								L,
								lengths.count(L)
								);
						}
					
					int count_out=0;
					for(String dna:dnaIndexes.keySetDecreasing())
						{
						if(++count_out>owner().COUNT_INDEX) break;
						tsv(owner().wDNAIndexes,this.hash,dna,dnaIndexes.count(dna));
						}
					}
				
			}
		}

	private int COUNT_INDEX=0;
	private File OUT=null;
	
	@Override
	public void printOptions(PrintStream out) {
		out.println(" -h help. This Screen");
		out.println(" -v print version and exits");
		out.println(" -L (log-level , a "+Level.class.getName()+"). Optional");
		out.println(" -X (int) maximum number of DNA indexes to print. memory consuming if not 0. Optional");
		out.println(" -o (filename out) Directory or ZIP. Required.");
		out.println(" -T (int) number of threads current:"+this.nThreads);
		}
	
	@Override
	public int doWork(String args[])
		{
		GetOpt getopt=new GetOpt();
		int c;
		while((c=getopt.getopt(args, "hvL:o:X:T:"))!=-1)
			{
			switch(c)
				{
				case 'h': printUsage();return 0;
				case 'v': System.out.println(getVersion());return 0;
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
				case 'T':
					{
					this.nThreads=Math.max(1, Integer.parseInt(getopt.getOptArg()));
					break;
					}
				case ':':
					{
					System.err.println("Missing argument. -"+(char)getopt.getOptOpt());
					return -1;
					}
				default:
					{
					System.err.println("Unknown option -"+(char)getopt.getOptOpt());
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
			this.wsqlite = archiveFactory.openWriter("sqlite3.sql");
			
			info("reading from stdin");
			BufferedReader in=new BufferedReader(new InputStreamReader(System.in));
			List<Analyzer> pool=new ArrayList<Analyzer>(this.nThreads);
			for(;;)
				{
				String line=in.readLine();
				if(line==null || pool.size()==this.nThreads)
					{
					for(Analyzer analyzer:pool)
						{
						analyzer.start();
						}
					for(Analyzer analyzer:pool)
						{
						analyzer.join();
						}
					pool.clear();
					if(line==null) break;
					}
				pool.add(new Analyzer(new File(line)));
				}
			
			in.close();
			
		
			this.wsqlite.println("create table if not exists wnames ( path TEXT, directory PATH, filename TEXT, md5 TEST,sample TEXT, dnaIndex TEXT, lane INT, side TEXT, split INT, fileSize INT );");
			this.wsqlite.println("create table if not exists wcount ( md5 TEXT, nReads INT, count_read_fails_filter INT, count_read_doesnt_fail_filter INT );");
			this.wsqlite.println("create table if not exists wquals ( md5 TEXT, qual FLOAT );");
			this.wsqlite.println("create table if not exists whistquals ( md5 TEXT, qual FLOAT, nqual INT );");
			this.wsqlite.println("create table if not exists wqualperpos ( md5 TEXT, position INT, qual FLOAT, nbases INT );");
			this.wsqlite.println("create table if not exists wbases ( md5 TEXT, position INT, A INT, T INT, G INT, C INT, N INT );");
			this.wsqlite.println("create table if not exists wlength ( md5 TEXT, readLen INT, nReads INT );");
			this.wsqlite.println("create table if not exists wDNAIndexes ( md5 TEXT, dnaIndex TEXT, nReads INT );");
			this.wsqlite.println(".separator '\t'");
			
			this.wsqlite.println(".import  names.tsv  wnames\n");
			this.wsqlite.println(".import  counts.tsv  wcount\n");
			this.wsqlite.println(".import  quals.tsv  wquals\n");
			this.wsqlite.println(".import  histquals.tsv  whistquals\n");
			this.wsqlite.println(".import  histpos2qual.tsv  wqualperpos\n");
			this.wsqlite.println(".import  bases.tsv  wbases\n");
			this.wsqlite.println(".import  lengths.tsv  wlength\n");
			this.wsqlite.println(".import  indexes.tsv  wDNAIndexes\n");
			
			
			for(PrintWriter pw: new PrintWriter[]{
					this.wnames,
					this.wcount,
					this.wquals,
					this.wbadfastq,
					this.whistquals,
					this.wqualperpos,
					this.wbases,
					this.wlength,
					this.wDNAIndexes,
					this.wsqlite
					})
				{
				pw.flush();
				pw.close();
				}
			

			
			if(archiveFactory!=null)
				{
				this.archiveFactory.close();
				}
			info("Done.");
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
