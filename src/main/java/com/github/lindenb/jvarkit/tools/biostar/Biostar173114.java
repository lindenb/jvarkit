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

package com.github.lindenb.jvarkit.tools.biostar;
import java.util.Collection;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;


public class Biostar173114 extends AbstractBiostar173114
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(Biostar173114.class);

	
	
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		SAMRecordIterator iter=null;
		try
			{
			BlockCompressedOutputStream.setDefaultCompressionLevel(9);
			sfr = super.openSamReader(inputName);
			final SAMFileHeader header0 = sfr.getFileHeader();
			final SAMFileHeader header = new SAMFileHeader();
			header.setSequenceDictionary(header0.getSequenceDictionary());
			header.setGroupOrder(header0.getGroupOrder());
			header.setSortOrder(header0.getSortOrder());
			sfw = super.openSAMFileWriter(header, true);
			iter=sfr.iterator();
			final SAMRecordFactory samRecordFactory = new DefaultSAMRecordFactory();
			long nReads=0;
			while(iter.hasNext())
				{
				final SAMRecord record= iter.next();
				
				
				if(record.getReadUnmappedFlag()) continue;
				if(record.getReadFailsVendorQualityCheckFlag()) continue;
				if(record.isSecondaryOrSupplementary()) continue;
				if(record.getDuplicateReadFlag()) continue;
				if(record.getMappingQuality()< super.minMapQ) continue;
				
				final SAMRecord rec2=samRecordFactory.createSAMRecord(header);
				 rec2.setReadUnmappedFlag(false);
				 rec2.setReadNegativeStrandFlag(record.getReadNegativeStrandFlag());;
				 rec2.setReadName("R"+ Long.toHexString(nReads++));
				 rec2.setReferenceIndex(record.getReferenceIndex());
				 rec2.setAlignmentStart(record.getAlignmentStart());
				 if(record.getCigar()!=null) {
				 final Cigar cigar = record.getCigar();
				 final java.util.List<CigarElement> cl = new java.util.ArrayList<>(cigar.numCigarElements());
				 for(int i=0;i<cigar.numCigarElements();++i) {
				  final CigarElement ce=cigar.getCigarElement(i);
				  if(ce.getOperator()==CigarOperator.S || ce.getOperator()==CigarOperator.H)
				    {
				    continue;
				    }
				  cl.add(ce);
				  }
				 
				 rec2.setCigar(new Cigar(cl));
				 }
				 
				 rec2.setMappingQuality(record.getMappingQuality());
				 rec2.setReadBases(super.rmSeq?SAMRecord.NULL_SEQUENCE:record.getReadBases());
				 rec2.setBaseQualities(SAMRecord.NULL_QUALS);
				sfw.addAlignment(rec2);
				}
			LOG.info("done");
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(sfr);
			CloserUtil.close(sfw);
			}
		}
	
	public static void main(String[] args)throws Exception
		{
		new Biostar173114().instanceMain(args);
		}
	}
