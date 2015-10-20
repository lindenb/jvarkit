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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloserUtil;

import java.util.Collection;
import java.util.List;


import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class SamRetrieveSeqAndQual extends AbstractSamRetrieveSeqAndQual
	{

	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(SamRetrieveSeqAndQual.class);

	 @Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  class MyCommand extends AbstractSamRetrieveSeqAndQual.AbstractSamRetrieveSeqAndQualCommand
	 	{

	
	private String normalizeFastqName(String s)
		{
		int w=s.indexOf(' ');
		if(w!=-1) s=s.substring(0,w);
		if(s.endsWith("/1") || s.endsWith("/2")) s=s.substring(0,s.length()-2);
		return s;
		}
	
	@SuppressWarnings("resource")
	@Override
		protected Collection<Throwable> call(String inputName) throws Exception {
			FastqReader[] fastqReaders=null;
			SamReader samReader=null;
			SAMFileWriter samWriter=null;
			SAMRecordIterator iter=null;
			try
				{
				
				if(super.fastqFin==null)
					{
					return wrapException("undefined fastq file");
					}
				else 
					{
					LOG.info("opening "+fastqFin);
					FastqReader r1= openFastqFileReader(fastqFin);
					if(fastqRin==null)
						{
						fastqReaders=new FastqReader[]{r1};
						}
					else
						{
						LOG.info("opening "+fastqRin);
						FastqReader r2 = openFastqFileReader(fastqRin);
						fastqReaders=new FastqReader[]{r1,r2};
						}
					}
	
				
				
				
				
				samReader  = openSamReader(inputName);
					
				
				SAMFileHeader.SortOrder sortOrder = samReader.getFileHeader().getSortOrder();
				if(sortOrder==null)
					{
					LOG.warn("undefined sort order read are in the sam order");
					}
				else if(sortOrder.equals(SAMFileHeader.SortOrder.coordinate))
					{
					return wrapException("Bad Sort Order. Sort this input on read name");
					}
				
				SAMFileHeader header= samReader.getFileHeader().clone();
				samWriter = openSAMFileWriter(header, true);
					
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
							return wrapException("Not paired but number of fastq!=2");
							}
						fastq_index = (rec.getFirstOfPairFlag()?0:1);
						}
					else
						{
						if(fastqReaders.length!=1)
							{
							return wrapException("Not paired but number of fastq!=1");
							}
						fastq_index=0;
						}
					
					if(sortOrder==SAMFileHeader.SortOrder.queryname)
						{
						while(currFastq[fastq_index]==null || normalizeFastqName(currFastq[fastq_index].getReadHeader()).compareTo(readName)<0)
							{
							if(! fastqReaders[ fastq_index ].hasNext())
								{
								return wrapException("Read Missing for "+readName);
								}
							currFastq[fastq_index] = fastqReaders[ fastq_index ].next();
							if(normalizeFastqName(currFastq[fastq_index].getReadHeader()).compareTo(readName)>0)
								{
								return wrapException("Read Missing for "+readName);
								}
							}
						}
					else
						{
						if(! fastqReaders[ fastq_index ].hasNext())
							{
							return wrapException("Read Missing for "+readName);
							}
						currFastq[fastq_index] = fastqReaders[ fastq_index ].next();
						}
					if(normalizeFastqName(currFastq[fastq_index].getReadHeader()).compareTo(readName)!=0)
						{
						return wrapException("Read Missing/Error for "+readName+" current:" + currFastq[fastq_index].getReadHeader());
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
				return RETURN_OK;
				}
			catch (Exception err)
				{
				return wrapException(err);
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
	 	}

	public static void main(String[] args)
		{
		new SamRetrieveSeqAndQual().instanceMainWithExit(args);
		}

	}
