/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


History:
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloserUtil;

import java.io.File;
import java.util.List;


import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class SamRetrieveSeqAndQual extends AbstractCommandLineProgram
	{

	@Override
	public String getProgramDescription() {
		return "I have a query-sorted BAM file without read/qual sequences and a FASTQ file with the read/qual sequences. Is there a tool to add seq to BAM?  for @sjackman https://twitter.com/sjackman/status/575368165531611136";
		}			
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/SamRetrieveSeqAndQual";
		}
	
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -o (file.out) optional. default: stdout");
		out.println(" -F (fastq / fastqF) required.");
		out.println(" -R (fastqR ) required.");
		super.printOptions(out);
		}
	
	private String normalizeFastqName(String s)
		{
		int w=s.indexOf(' ');
		if(w!=-1) s=s.substring(0,w);
		if(s.endsWith("/1") || s.endsWith("/2")) s=s.substring(0,s.length()-2);
		return s;
		}
	
	@Override
	public int doWork(String[] args)
		{
		
		File fastqFin=null;
		File fastqRin=null;
		File bamOut=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:F:R:"))!=-1)
			{
			switch(c)
				{
				case 'o': bamOut=new File(opt.getOptArg()); break;
				case 'F': fastqFin=new File(opt.getOptArg()); break;
				case 'R': fastqRin=new File(opt.getOptArg()); break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
			
		FastqReader[] fastqReaders=null;
		SamReader samReader=null;
		SAMFileWriter samWriter=null;
		SAMRecordIterator iter=null;
		try
			{
			
			if(fastqFin==null)
				{
				error("undefined fastq file");
				return -1;
				}
			else 
				{
				info("opening "+fastqFin);
				FastqReader r1=new FastqReader(fastqFin);
				if(fastqRin==null)
					{
					fastqReaders=new FastqReader[]{r1};
					}
				else
					{
					info("opening "+fastqRin);
					FastqReader r2=new FastqReader(fastqRin);
					fastqReaders=new FastqReader[]{r1,r2};
					}
				}

			
			
			SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			
			int optind = opt.getOptInd();
			if(optind==args.length)
				{
				info("Reading BAM from stdin");
				samReader  = srf.open(SamInputResource.of(System.in));
				}
			else if(optind+1==args.length)
				{
				info("Reading BAM from "+args[optind]);
				samReader  = srf.open(new File(args[optind]));
				}
			else
				{
				error("illegal.number.of.arguments");
				return -1;
				}
			
			SAMFileHeader.SortOrder sortOrder = samReader.getFileHeader().getSortOrder();
			if(sortOrder==null)
				{
				warning("undefined sort order read are in the sam order");
				}
			else if(sortOrder.equals(SAMFileHeader.SortOrder.coordinate))
				{
				error("Bad Sort Order. Sort this input on read name");
				return -1;
				}
			
			SAMFileHeader header= samReader.getFileHeader().clone();
			
			SAMProgramRecord prg=header.createProgramRecord();
			prg.setCommandLine(this.getProgramCommandLine());
			prg.setProgramName(this.getProgramName());
			prg.setProgramVersion(this.getVersion());
			
			SAMFileWriterFactory sfwf=new SAMFileWriterFactory();
			sfwf.setCreateIndex(false);
			sfwf.setCreateMd5File(false);
			if(bamOut!=null)
				{
				samWriter = sfwf.makeSAMOrBAMWriter(header, true, bamOut);
				}
			else
				{
				samWriter =sfwf.makeSAMWriter(header, true, System.out);
				}
			iter=samReader.iterator();
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			FastqRecord currFastq[]=new FastqRecord[]{null,null};
			while(iter.hasNext())
				{
				SAMRecord rec= progress.watch(iter.next());
				String readName = rec.getReadName();
				int fastq_index = 0;
				if(rec.getReadPairedFlag())
					{
					if(fastqReaders.length!=2)
						{
						error("Not paired but number of fastq!=2");
						return  -1;
						}
					fastq_index = (rec.getFirstOfPairFlag()?0:1);
					}
				else
					{
					if(fastqReaders.length!=1)
						{
						error("Not paired but number of fastq!=1");
						return  -1;
						}
					fastq_index=0;
					}
				
				if(sortOrder==SAMFileHeader.SortOrder.queryname)
					{
					while(currFastq[fastq_index]==null || normalizeFastqName(currFastq[fastq_index].getReadHeader()).compareTo(readName)<0)
						{
						if(! fastqReaders[ fastq_index ].hasNext())
							{
							error("Read Missing for "+readName);
							return -1;
							}
						currFastq[fastq_index] = fastqReaders[ fastq_index ].next();
						if(normalizeFastqName(currFastq[fastq_index].getReadHeader()).compareTo(readName)>0)
							{
							error("Read Missing for "+readName);
							return -1;
							}
						}
					}
				else
					{
					if(! fastqReaders[ fastq_index ].hasNext())
						{
						error("Read Missing for "+readName);
						return -1;
						}
					currFastq[fastq_index] = fastqReaders[ fastq_index ].next();
					}
				if(normalizeFastqName(currFastq[fastq_index].getReadHeader()).compareTo(readName)!=0)
					{
					error("Read Missing/Error for "+readName+" current:" + currFastq[fastq_index].getReadHeader());
					return -1;
					}
				String fastqBases = currFastq[fastq_index].getReadString();
				String fastqQuals = currFastq[fastq_index].getBaseQualityString();
				/* handle orientation */
				if(!rec.getReadUnmappedFlag() && rec.getReadNegativeStrandFlag())
					{
					fastqBases = AcidNucleics.reverseComplement(fastqBases);
					StringBuilder sb=new StringBuilder(fastqQuals.length());
					for(int i=fastqQuals.length()-1;i>=0;--i)
						sb.append(fastqQuals.charAt(i));
					fastqQuals =  sb.toString();
					}
				/* remove hard clip */
				Cigar cigar=rec.getCigar();
				if(cigar!=null)
					{
					List<CigarElement> ceList=cigar.getCigarElements();
					if(!ceList.isEmpty())
						{
						CigarElement ce = ceList.get(ceList.size()-1);
						if(ce.getOperator()==CigarOperator.HARD_CLIP)
							{
							fastqBases = fastqBases.substring(0, fastqBases.length()-ce.getLength());
							fastqQuals = fastqQuals.substring(0, fastqQuals.length()-ce.getLength());
							}
						ce = ceList.get(0);
						if(ce.getOperator()==CigarOperator.HARD_CLIP)
							{
							fastqBases = fastqBases.substring(ce.getLength());
							fastqQuals = fastqQuals.substring(ce.getLength());
							}
						}
					}
				rec.setBaseQualityString(fastqQuals);
				rec.setReadString(fastqBases);
				samWriter.addAlignment(rec);
				}
			progress.finish();
			return 0;
			}
		catch (Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			if(fastqReaders!=null)
				for(FastqReader r:fastqReaders)
					CloserUtil.close(r);
			CloserUtil.close(iter);
			CloserUtil.close(samReader);
			CloserUtil.close(samWriter);
			}
		}


	public static void main(String[] args)
		{
		new SamRetrieveSeqAndQual().instanceMainWithExit(args);
		}

	}
