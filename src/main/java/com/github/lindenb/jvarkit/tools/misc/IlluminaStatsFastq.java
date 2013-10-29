package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMUtils;

import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.illumina.FastQName;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;

public class IlluminaStatsFastq
	extends AbstractCommandLineProgram
	{
	private static final Log LOG=Log.getInstance(IlluminaStatsFastq.class);
	
	
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"  Count FASTQs in Illumina Result. Generate a script for sqlite3 ";
	@Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Directories to process",optional=false)
	public File IN=null;
	
	@Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output Name (Directory or .zip file)",optional=false)
	public File OUT=null;
	
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
    
	
	
	private static void tsv(PrintWriter out,Object...array)
		{
		for(int i=0;i< array.length;++i)
			{
			if(i>0) out.print('\t');
			out.print(array[i]);
			}
		out.println();
		}
	
	private void recursive(File f) throws IOException
		{
		if(f==null) return;
		LOG.info(f);
		if(f.isDirectory())
			{
			
			File children[]=f.listFiles();
			if(children==null) return;
			for(File c:children)
				{
				recursive(c);
				}
			} 
		else if(f.getName().endsWith(".fastq.gz") && f.isFile())
			{
			final int QUALITY_STEP=5;
			
			FastQName fq=FastQName.parse(f);
			
			Counter<Integer> qualityHistogram=new Counter<Integer>();
			Counter<Integer> pos2quality=new Counter<Integer>();
			List<Bases> pos2bases=new ArrayList<Bases>(300);
			Counter<Integer> lengths=new Counter<Integer>();
			Counter<Integer> pos2count=new Counter<Integer>();

			long nReads=0L;
			double sum_qualities=0L;
			long count_bases=0L;
			long count_read_fails_filter=0L;
			long count_read_doesnt_fail_filter=0L;
			FastqReader r=null;
			try
				{
				r=new FastqReader(f);
				}
			catch(Exception error)
				{
				error.printStackTrace();
				tsv(wbadfastq,f.getPath(),error.getMessage());
				return;
				}	
			
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
			r.close();
			
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
			
			
			}
		}
	
	
	@Override
	protected int doWork()
		{
		try {
			
			if(!IN.exists())
				{
				LOG.error("Input "+IN+" doesnt exists.");
				return -1;
				}
			
			archiveFactory=ArchiveFactory.open(OUT);
			this.wnames = archiveFactory.openWriter("names.tsv");
			this.wcount = archiveFactory.openWriter("counts.tsv");
			this.wquals = archiveFactory.openWriter("quals.tsv");
			this.wbadfastq = archiveFactory.openWriter("notfastq.tsv");
			this.whistquals = archiveFactory.openWriter("histquals.tsv");
			this.wqualperpos = archiveFactory.openWriter("histpos2qual.tsv");
			this.wbases = archiveFactory.openWriter("bases.tsv");
			this.wlength = archiveFactory.openWriter("lengths.tsv");
			
			recursive(IN);
			
			for(PrintWriter pw: new PrintWriter[]{
					this.wnames,
					this.wcount,
					this.wquals,
					this.wbadfastq,
					this.whistquals,
					this.wqualperpos,
					this.wbases,
					this.wlength})
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
			LOG.error(e, "Error");
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
