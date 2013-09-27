package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
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
			if(sql)
				{
				System.out.println(
						"insert into FASTQ(file,sample,seqindex,lane,side,split,size) values ("+
						quote(fq.getFile().getPath())+","+
						quote(fq.isUndetermined()?"Undetermined":fq.getSample())+","+
						quote(fq.getSeqIndex())+","+
						fq.getLane()+","+
						fq.getSide()+","+
						fq.getSplit()+","+
						fq.getFile().length()+
						");");
				//+ile\tsample\tseqindex\tlane\tside\tsplit\tsize
				}
			else
				{
				System.out.println(
					fq.getFile().getPath()+"\t"+
					(fq.isUndetermined()?"Undetermined":fq.getSample())+"\t"+
					fq.getLane()+"\t"+
					fq.getSide()+"\t"+
					fq.getSplit()+"\t"+
					fq.getFile().length()
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
				System.out.println("create table  if not exists FASTQ(file TEXT,sample text,seqindex text,lane int,side int,split int,size int);");
				System.out.println("begin transaction;");
				}
			else
				{
				System.out.println("#file\tsample\tseqindex\tlane\tside\tsplit\tsize");
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
