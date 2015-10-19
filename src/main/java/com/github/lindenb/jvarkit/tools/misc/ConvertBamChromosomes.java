/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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




*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.SAMFileHeader.SortOrder;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;


import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class ConvertBamChromosomes
	extends AbstractConvertBamChromosomes
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(ConvertBamChromosomes.class);

	 @Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  class MyCommand extends AbstractConvertBamChromosomes.AbstractConvertBamChromosomesCommand
	 	{
		private Map<String,String> customMapping=new HashMap<String,String>();
		private Set<String> unmappedChromosomes=new HashSet<String>();
		
		
		private String convertName(String chrom)throws IOException
			{
			if(chrom==null) throw new NullPointerException();
			String newname=customMapping.get(chrom);
			if(newname==null)
				{
				if(!unmappedChromosomes.contains(chrom))
					{
					LOG.warn("unmapped chromosome "+chrom);
					unmappedChromosomes.add(chrom);
					}
				if(ignore_if_no_mapping) return null;
				
				if(use_original_chrom_name_if_no_mapping)
					{	
					return chrom;
					}
				throw new IOException("No mapping found to convert name of chromosome \""+chrom+"\"");
				}
			return newname;
			}
		
		
		
		@Override
		protected Collection<Throwable> call(String inputName) throws Exception {
			if(super.mappingFile==null)
				{
				return wrapException("undefined mapping file");
				}
					
			BufferedReader in=null;
			try
				{
				customMapping = super.loadCustomChromosomeMapping(super.mappingFile);
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(in);
				}
					
			SamReader sfr=null;
			SAMFileWriter sfw=null;
			try
				{
				sfr  = openSamReader(inputName);
				SAMFileHeader header1=sfr.getFileHeader();
				if(header1==null)
					{
					return wrapException("File header missing");
					}
				
				SAMFileHeader header2=header1.clone();
				
				//create new sequence dict
				final SAMSequenceDictionary dict1=header1.getSequenceDictionary();
				if(dict1==null)
					{
					return wrapException("Sequence dict missing");
					}
				List<SAMSequenceRecord> ssrs=new ArrayList<SAMSequenceRecord>(dict1.size());
				for(int i=0;i< dict1.size();++i)
					{
					SAMSequenceRecord ssr=dict1.getSequence(i);
					String newName=convertName(ssr.getSequenceName());
					if(newName==null)
						{
						//skip unknown chromosomes
						continue;
						}
					ssr=new SAMSequenceRecord(newName, ssr.getSequenceLength());
					ssrs.add(ssr);
					}
				header2.setSequenceDictionary(new SAMSequenceDictionary(ssrs));
				
				SAMSequenceDictionary dict2=new SAMSequenceDictionary(ssrs);
				header2.setSequenceDictionary(dict2);
				SAMProgramRecord prog=header2.createProgramRecord();
				prog.setCommandLine(this.getProgramCommandLine());
				prog.setProgramName(getName());
				prog.setProgramVersion(getVersion());
				
				
				
				
				SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict1);
	
				boolean presorted=(header1.getSortOrder()!=null && (header1.getSortOrder()==SortOrder.coordinate || header1.getSortOrder()==SortOrder.queryname));
				sfw = openSAMFileWriter(header2, presorted);
				
				
				
				long num_ignored=0L;
				SAMRecordIterator iter=sfr.iterator();
				while(iter.hasNext())
					{
					SAMRecord rec1=iter.next();
					progress.watch(rec1);
					String newName1=null;
					String newName2=null;
					if(!SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getReferenceName()))
						{
						newName1=convertName(rec1.getReferenceName());
						}
					if(rec1.getReadPairedFlag() && !SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getMateReferenceName()))
						{
						newName2=convertName(rec1.getMateReferenceName());
						}
					rec1.setHeader(header2);
	
					if(!SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getReferenceName()))
						{
						if(newName1==null)
							{
							++num_ignored;
							continue;
							}
						rec1.setReferenceName(newName1);
						}
					if(rec1.getReadPairedFlag() && !SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getMateReferenceName()))
						{
						if(newName2==null)
							{
							++num_ignored;
							continue;
							}
						rec1.setMateReferenceName(newName2);
						}
					sfw.addAlignment(rec1);
					}
				if(!unmappedChromosomes.isEmpty())
					{
					LOG.warn("Unmapped chromosomes: "+unmappedChromosomes);
					}
				LOG.warn("num ignored read:"+num_ignored);
				return RETURN_OK;
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(sfr);
				CloserUtil.close(sfw);
				}
			}
	
	 	}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new ConvertBamChromosomes().instanceMainWithExit(args);
		}
	}
