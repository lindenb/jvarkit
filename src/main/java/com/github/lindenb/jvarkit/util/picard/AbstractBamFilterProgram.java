package com.github.lindenb.jvarkit.util.picard;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;

import com.github.lindenb.jvarkit.util.cli.GetOpt;


import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.util.CloserUtil;


public abstract class AbstractBamFilterProgram
	extends com.github.lindenb.jvarkit.util.AbstractCommandLineProgram
		{
		protected File bamFileOut=null;
		protected ValidationStringency validationStringency=ValidationStringency.LENIENT;
		
		@Override
		public String getProgramDescription() {
			return "Filter Sam";
			}
		
	
		protected boolean create_md5=false;
		protected boolean create_index=false;
		protected int max_records_in_ram=50000;

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
			String vers=getVersion().trim();
			if(!vers.isEmpty()) spr.setProgramVersion(vers);
			//h2.addProgramRecord(spr);
			return h2;
			}
		
		
		@Override
		public void printOptions(PrintStream out)
			{
			out.println(" -o (filename) out. Optional");
			out.println(" -V (Validation stringency ) Optional. Current is :"+validationStringency);
			out.println(" -c (int) max records in ram. Optional. Current is :"+max_records_in_ram);
			super.printOptions(out);
			}
		
		@Override
		protected String getGetOptDefault()
			{
			return super.getGetOptDefault()+"V:o:c:";
			}
		
		@Override
		protected GetOptStatus handleOtherOptions(int c, GetOpt opt, String[] args)
			{
			switch(c)
				{
				case 'c':
					{
					try
						{
						this.max_records_in_ram=Integer.parseInt(opt.getOptArg());
						}
					catch(NumberFormatException err)
						{
						error("Bad value for max_records_in_ram");
						return GetOptStatus.EXIT_FAILURE;
						}
					return GetOptStatus.OK;
					}
				case 'o':
					{
					this.bamFileOut=new File(opt.getOptArg());
					return GetOptStatus.OK;
					}
				case 'V':
					{
					try
						{
						this.validationStringency=ValidationStringency.valueOf(opt.getOptArg());
						}
					catch(Exception err)
						{
						error("Bad validation stringency: expected one of: "+Arrays.asList(ValidationStringency.values()));
						return GetOptStatus.EXIT_FAILURE;
						}
					return GetOptStatus.OK;
					}
				default:return super.handleOtherOptions(c, opt, args);
				}
			
			}
		
		protected boolean isOutputPresorted(final SAMFileHeader header)
			{
			return	header!=null &&
					header.getSortOrder()!=null &&
					header.getSortOrder()!=SortOrder.unsorted
					;
			}
		
		protected int doWork(File fileIn,File fileOut)
			throws Exception
			{
			SAMFileReader sfr=null;
			SAMFileWriter sfw=null;
			SAMFileHeader header;
			try
				{

				
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
				
				sfr.setValidationStringency(this.validationStringency);
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
							isOutputPresorted(headerOut)
							,System.out);
					}
				else
					{
					info("Writing to "+fileOut);
					sfw=sfwf.makeSAMOrBAMWriter(
							headerOut, 
							isOutputPresorted(headerOut)
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
