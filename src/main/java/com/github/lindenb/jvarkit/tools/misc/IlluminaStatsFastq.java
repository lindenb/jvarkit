package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.util.Log;

import com.github.lindenb.jvarkit.util.illumina.FastQName;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;

public class IlluminaStatsFastq
	extends AbstractCommandLineProgram
	{
	private static final Log LOG=Log.getInstance(IlluminaStatsFastq.class);
	
	
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"  Count FASTQs in Illumina Result";
	@Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Directories to process",optional=false)
	public File IN=null;

	@Option(shortName="SQLITE", doc="sqlite output",optional=false)
	public boolean sql=false;
	
    @Option(shortName="CAT",doc="category name",optional=true)
    public String CATEGORY="undefined";

	
   
    
	private static String quote(String s)
		{
		return "'"+s+"'";
		}
	
	private void recursive(File f) throws IOException
		{
		if(f==null) return;
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
			FastQName fq=FastQName.parse(f);
			if(!fq.isValid())
				{
				LOG.info("invalid file");
				return;
				}
			
			Map<Integer,Integer> lengths2count=new TreeMap<Integer, Integer>();
			List<long[]> pos2bases=new ArrayList<long[]>(100);
			List<Map<Integer,Long>> pos2qual=new ArrayList<Map<Integer,Long>>(100);
			Map<String,Integer> indexseq2count=new HashMap<String, Integer>();
			
			long nReads=0L;
			
			
			FastqReader fqr=null;
			try
				{
				fqr=new FastqReader(f);
				while(fqr.hasNext())
					{
					FastqRecord r=fqr.next();
					++nReads;
					int readLen=r.getReadString().length();
					Integer countLength=lengths2count.get(readLen);
					lengths2count.put(readLen,
							countLength==null?1:countLength+1
							);
					
					while(pos2bases.size()< readLen)
						{
						pos2bases.add(new long[]{0L,0L,0L,0L,0L});
						pos2qual.add(new HashMap<Integer,Long>(128));
						}
					
					
					for(int b=0;b< readLen ;++b)
						{
						int countBaseIndex=-1;
						switch(r.getReadString().charAt(b))
							{
							case 'A': case 'a':countBaseIndex=0;break;
							case 'T': case 't':countBaseIndex=1;break;
							case 'G': case 'g':countBaseIndex=2;break;
							case 'C': case 'c':countBaseIndex=3;break;
							default: countBaseIndex=4;break;
							}
						pos2bases.get(b)[countBaseIndex]++;
						}
					
					for(int b=0;b< readLen ;++b)
						{
						int qual=(int)r.getBaseQualityString().charAt(b) - 33;
						Map<Integer,Long> qual2count=pos2qual.get(b);
						Long qualc=qual2count.get(qual);
						qual2count.put(qual,
								qualc==null?1:qualc+1
								);
						}
					
					
					IlluminaReadName fqName=IlluminaReadName.parse1_4(r.getReadHeader());
					String indexSeq=fqName.getIndex();
					Integer countIndex=indexseq2count.get(indexSeq);
					indexseq2count.put(indexSeq,
							countIndex==null?1:countIndex+1
							);
					}
				

				
				
				}
			catch(Exception err)
				{
				nReads=-1L;
				err.printStackTrace();
				}
			finally
				{
				if(fqr!=null) fqr.close();
				}
			
			
			if(sql)
				{
				System.out.println(
						"insert into FASTQ(category,file,sample,seqindex,lane,side,split,size,countReads) values ("+
						quote(this.CATEGORY)+","+		
						quote(fq.getFile().getPath())+","+
						quote(fq.isUndetermined()?"Undetermined":fq.getSample())+","+
						quote(fq.getSeqIndex())+","+
						fq.getLane()+","+
						fq.getSide()+","+
						fq.getSplit()+","+
						fq.getFile().length()+","+
						nReads+
						");");

				
				for(String seqindex:indexseq2count.keySet())
					{
					System.out.println(
							"insert into SEQINDEX(file,seqindex,countindex) values ("+
							quote(fq.getFile().getPath())+","+
							quote(seqindex)+","+
							indexseq2count.get(seqindex)+
							");");
					}
				
				for(int b=0;b< pos2bases.size();++b)
					{
					long counts[]=pos2bases.get(b);
					System.out.println(
							"insert into BASES(file,position,A,T,G,C,N) values ("+
							quote(fq.getFile().getPath())+","+
							(b+1)+","+
							counts[0]+","+counts[1]+","+counts[2]+","+counts[3]+","+counts[4]+
							");");
					}
				for(int b=0;b< pos2qual.size();++b)
					{
					Map<Integer,Long> qual2count=pos2qual.get(b);
					long total=0L;
					long count=0L;
					for(Integer q:qual2count.keySet())
						{	
						long n_qual=qual2count.get(q);
						total+=q*n_qual;
						count+=n_qual;
						}
					double mean=(total)/(double)count;
					System.out.println(
							"insert into QUALITY(file,position,qual) values ("+
							quote(this.CATEGORY)+","+		
							(b+1)+","+
							mean+
							");");
					}
				}
			else
				{
				System.out.println(
					this.CATEGORY+"\t"+
					fq.getFile().getPath()+"\t"+
					(fq.isUndetermined()?"Undetermined":fq.getSample())+"\t"+
					fq.getLane()+"\t"+
					fq.getSide()+"\t"+
					fq.getSplit()+"\t"+
					fq.getFile().length()+"\t"+
					nReads
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
			if(!IN.isDirectory())
				{
				LOG.error("Input "+IN+" is not a directory.");
				return -1;
				}
			if(sql)
				{
				System.out.println("create table  if not exists FASTQ(category TEXT,file TEXT,sample text,seqindex text,lane int,side int,split int,size int,countReads int);");
				System.out.println("create table  if not exists SEQINDEX(file TEXT,seqindex text,int countindex);");
				System.out.println("create table  if not exists BASES(file TEXT,position int,A int,T int,G int,C int,N int);");
				System.out.println("create table  if not exists QUALITY(file TEXT,position int,qual float);");
				System.out.println("begin transaction;");
				}
			else
				{
				System.out.println("#category\tfile\tsample\tseqindex\tlane\tside\tsplit\tsize\tcount-reads");
				}
			recursive(IN);
			if(sql)
				{
				System.out.println("commit;");
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
