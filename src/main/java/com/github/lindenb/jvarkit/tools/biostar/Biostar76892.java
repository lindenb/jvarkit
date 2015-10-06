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
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;


public class Biostar76892 extends AbstractBiostar76892
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Biostar76892.class);

	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	
	static private class MyCommand extends AbstractBiostar76892.AbstractBiostar76892Command
		{
		@Override
		protected Collection<Throwable> call(String inputName) throws Exception {
			SamReader sfr=null;
			SAMFileWriter sfw=null;
			try
				{
				sfr= super.openSamReader(inputName);
				SAMFileHeader header=sfr.getFileHeader();
				if(header.getSortOrder()!=SortOrder.coordinate)
					{
					return wrapException("input must be sorted on coordinate ? got: "+header.getSortOrder() );
					}
				
				header= header.clone();
				header.addComment(getProgramCommandLine());
				sfw = super.openSAMFileWriter(header, true);
				SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header);
				List<SAMRecord> buffer=new ArrayList<SAMRecord>();
				SAMRecordIterator iter=sfr.iterator();
				for(;;)
					{
					SAMRecord rec=null;
					//get next record
					if(iter.hasNext())
						{
						rec=progress.watch(iter.next());
						if(!rec.getReadPairedFlag() ||
							rec.getReadUnmappedFlag() ||
							rec.getMateUnmappedFlag() ||
							rec.getProperPairFlag() ||
							rec.getReferenceIndex()!=rec.getMateReferenceIndex() ||
							rec.getReadNegativeStrandFlag()==!rec.getMateNegativeStrandFlag()
							)
							{
							if(!onlySaveFixed) sfw.addAlignment(rec);
							continue;
							}
						}
	
					if(rec!=null)
						{
						int i=0;
						//cleanup buffer
						int mate_index=-1;
						while(i<buffer.size())
							{
							SAMRecord prev=buffer.get(i);
							if(prev.getReferenceIndex()!=rec.getReferenceIndex() ||
								prev.getAlignmentEnd() + distance < rec.getAlignmentStart())
								{
								if(!onlySaveFixed) sfw.addAlignment(prev);
								buffer.remove(i);
								}
							else if(prev.getReadName().equals(rec.getReadName()) && 
									(	(prev.getFirstOfPairFlag() && rec.getSecondOfPairFlag()) ||
										(rec.getFirstOfPairFlag() && prev.getSecondOfPairFlag()))
									)
								{
								mate_index=i;
								++i;
								}
							else
								{
								++i;
								}
							}
						
						if(mate_index==-1)
							{
							buffer.add(rec);
							}
						else
							{
							
							SAMRecord mate=buffer.get(mate_index);
							buffer.remove(mate_index);
							LOG.info("changing "+rec+" "+mate);
							
							if(mate.getReadNegativeStrandFlag())
								{
								mate.setReadNegativeStrandFlag(false);
								rec.setMateNegativeStrandFlag(mate.getReadNegativeStrandFlag());
								}
							else
								{
								rec.setReadNegativeStrandFlag(false);
								mate.setMateNegativeStrandFlag(rec.getReadNegativeStrandFlag());
								}
							if(!mate.getReadFailsVendorQualityCheckFlag() && !rec.getReadFailsVendorQualityCheckFlag())
								{
								mate.setProperPairFlag(true);
								rec.setProperPairFlag(true);
								}
							mate.setAttribute("rv",1);
							rec.setAttribute("rv", 1);
							sfw.addAlignment(mate);
							sfw.addAlignment(rec);
							}
						
						}
					else
						{
						for(SAMRecord r:buffer)
							{
							if(!onlySaveFixed)  sfw.addAlignment(r);
							}
						break;
						}
					}
				LOG.info("done");
				return Collections.emptyList();
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(sfw);
				CloserUtil.close(sfr);
				}
			}
		}
	public static void main(String[] args)throws Exception
		{
		new Biostar76892().instanceMain(args);
		}
	}
