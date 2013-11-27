package com.github.lindenb.jvarkit.util.picard;

import java.io.File;
import java.util.ArrayList;
import java.util.List;


import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.IOUtil;


public abstract class AbstractBamFilterProgram
	extends com.github.lindenb.jvarkit.util.AbstractCommandLineProgram
		{
		protected List<File> TMP_DIR=new ArrayList<File>();
	CommandLineProgram x;
		@Override
		public String getProgramDescription() {
			return "Filter Sam";
			}
		
	
		protected boolean create_md5=false;
		protected boolean create_index=false;
		protected int max_records_in_ram=5000;

		protected abstract int doWork(SAMFileReader in,SAMFileWriter out);
		
		
		protected boolean verifySortOrder(SortOrder sortOrder)
			{
			return true;
			}
		
		protected int doWork(File fileOut,int optind,String args[])
			{
			try
				{
				if(optind==args.length)
					{
					return doWork(null, fileOut);
					}
				else if(optind+1==args.length)
					{
					return doWork(new File(args[optind]), fileOut);
					}
				else
					{
					error("Illegal number of arguments.");
					return -1;
					}
				}
			catch(Exception err)
				{
				error(err);
				return -1;
				}
			}
		
		protected SortOrder createSamSortOrder(SortOrder inpuSortOrder)
			{
			if(inpuSortOrder==null) return SortOrder.unsorted;
			return inpuSortOrder;
			}
		
		protected SAMFileHeader createOuputSamFileHeader(SAMFileHeader inHeader)
			{
			SAMFileHeader h2= inHeader.clone();
			h2.setSortOrder(createSamSortOrder(inHeader.getSortOrder()));
			SAMProgramRecord spr=h2.createProgramRecord();
			spr.setProgramName(this.getProgramName());
			String cmd=getProgramCommandLine().replaceAll("[\t \n]+"," ").trim();
			if(!cmd.isEmpty()) spr.setCommandLine(cmd);
			h2.addProgramRecord(spr);
			return h2;
			}
		
		protected int doWork(File fileIn,File fileOut)
			throws Exception
			{
			SAMFileReader sfr=null;
			SAMFileWriter sfw=null;
			SAMFileHeader header;
			try
				{
		        if (this.TMP_DIR.isEmpty())
		        	{
		        	TMP_DIR.add(IOUtil.getDefaultTmpDir());
		        	}

				
				if( fileIn==null)
					{
					info("Reading stdin");
					sfr=new SAMFileReader(System.in);
					}
				else 
					{
					info("Reading "+fileIn);					
					sfr=new SAMFileReader(fileIn);
					}
				
				sfr.setValidationStringency(ValidationStringency.LENIENT);
				header=sfr.getFileHeader();
				if(!verifySortOrder(header.getSortOrder()))
					{
					return -1;
					}
				
				SAMFileWriterFactory sfwf=new SAMFileWriterFactory();
				sfwf.setCreateIndex(this.create_index);
				sfwf.setCreateMd5File(this.create_md5);
				sfwf.setMaxRecordsInRam(this.max_records_in_ram);
				
				SAMFileHeader headerOut=createOuputSamFileHeader(header);
				if(fileOut==null)
					{
					info("Writing to stdout");
					sfw=sfwf.makeSAMWriter(
							headerOut, 
							headerOut.getSortOrder()!=SortOrder.unsorted
							,System.out);
					}
				else
					{
					info("Writing to "+fileOut);
					sfw=sfwf.makeSAMOrBAMWriter(
							headerOut, 
							headerOut.getSortOrder()!=SortOrder.unsorted
							,fileOut);
					}
				
				return  doWork(sfr,sfw);
				}
			catch(Exception err)
				{
				error(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(sfr);
				CloserUtil.close(sfw);
				}
			}
	
		
}
