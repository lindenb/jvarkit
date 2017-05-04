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
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloserUtil;

import java.io.File;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

@Program(name="samretrieveseqandqual",description="I have a query-sorted BAM file without read/qual sequences and a FASTQ file with the read/qual sequences. Is there a tool to add seq to BAM?  for @sjackman https://twitter.com/sjackman/status/575368165531611136")
public class SamRetrieveSeqAndQual extends Launcher
	{
	private static final Logger LOG = Logger.build(SamRetrieveSeqAndQual.class).make();


	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File bamOut = null;
	@Parameter(names={"-F"},description=" (fastq / fastqF) required",required=true)
	private File fastqFin = null;
	@Parameter(names={"-R"},description=" (fastq / fastqR) required",required=true)
	private File fastqRin = null;

	@Parameter
	private WritingBamArgs writingBamArgs = new WritingBamArgs();
	
	private String normalizeFastqName(String s)
		{
		int w=s.indexOf(' ');
		if(w!=-1) s=s.substring(0,w);
		if(s.endsWith("/1") || s.endsWith("/2")) s=s.substring(0,s.length()-2);
		return s;
		}
	
	@Override
	public int doWork(List<String> args) {
		FastqReader[] fastqReaders=null;
		SamReader samReader=null;
		SAMFileWriter samWriter=null;
		SAMRecordIterator iter=null;
		try
			{
			
			if(fastqFin==null)
				{
				LOG.error("undefined fastq file");
				return -1;
				}
			else 
				{
				LOG.info("opening "+fastqFin);
				FastqReader r1=new FastqReader(fastqFin);
				if(fastqRin==null)
					{
					fastqReaders=new FastqReader[]{r1};
					}
				else
					{
					LOG.info("opening "+fastqRin);
					FastqReader r2=new FastqReader(fastqRin);
					fastqReaders=new FastqReader[]{r1,r2};
					}
				}

			samReader = super.openSamReader(oneFileOrNull(args));
			final SAMFileHeader.SortOrder sortOrder = samReader.getFileHeader().getSortOrder();
			if(sortOrder==null)
				{
				LOG.warning("undefined sort order read are in the sam order");
				}
			else if(!sortOrder.equals(SAMFileHeader.SortOrder.queryname))
				{
				LOG.error("Bad Sort Order. Sort this input on read name");
				return -1;
				}
			
			SAMFileHeader header= samReader.getFileHeader().clone();
			
			SAMProgramRecord prg=header.createProgramRecord();
			prg.setCommandLine(this.getProgramCommandLine());
			prg.setProgramName(this.getProgramName());
			prg.setProgramVersion(this.getVersion());
			
			samWriter  = writingBamArgs.openSAMFileWriter(bamOut, header, true);
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
						LOG.error("Not paired but number of fastq!=2");
						return  -1;
						}
					fastq_index = (rec.getFirstOfPairFlag()?0:1);
					}
				else
					{
					if(fastqReaders.length!=1)
						{
						LOG.error("Not paired but number of fastq!=1");
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
							LOG.error("Read Missing for "+readName);
							return -1;
							}
						currFastq[fastq_index] = fastqReaders[ fastq_index ].next();
						if(normalizeFastqName(currFastq[fastq_index].getReadHeader()).compareTo(readName)>0)
							{
							LOG.error("Read Missing for "+readName);
							return -1;
							}
						}
					}
				else
					{
					if(! fastqReaders[ fastq_index ].hasNext())
						{
						LOG.error("Read Missing for "+readName);
						return -1;
						}
					currFastq[fastq_index] = fastqReaders[ fastq_index ].next();
					}
				if(normalizeFastqName(currFastq[fastq_index].getReadHeader()).compareTo(readName)!=0)
					{
					LOG.error("Read Missing/Error for "+readName+" current:" + currFastq[fastq_index].getReadHeader());
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
			LOG.error(err);
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
