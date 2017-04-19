/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;

@Program(name="biostar76892",description="fix strand of two paired reads close but on the same strand. See http://www.biostars.org/p/76892/ ")
public class Biostar76892 extends Launcher
	{

	private static final Logger LOG = Logger.build(Biostar76892.class).make();
	
	
	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;
	
	
	@Parameter(names={"-osf","--osf"},description="only save pairs of reads which have been corrected by this program")
	private boolean onlySaveFixed = false;
	
	@Parameter(names={"-d","--maxc"},description="distance beween two reads")
	private int distance = 30 ;

	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();
    
	@Override
	public int doWork(final List<String> args) {
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		try
			{
			sfr= super.openSamReader(oneFileOrNull(args));
			sfw = writingBamArgs.openSAMFileWriter(this.outputFile,sfr.getFileHeader(), true);
			long nRecords=0;
			final List<SAMRecord> buffer=new ArrayList<SAMRecord>();
			SAMRecordIterator iter=sfr.iterator();
			for(;;)
				{
				SAMRecord rec=null;
				//get next record
				if(iter.hasNext())
					{
					rec=iter.next();
					++nRecords;
					if(nRecords%1000000==0) LOG.info("records: "+nRecords);
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
						final SAMRecord mate=buffer.get(mate_index);
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
					for(final SAMRecord r:buffer)
						{
						if(!onlySaveFixed)  sfw.addAlignment(r);
						}
					break;
					}
				}
			LOG.info("done");
			sfw.close();
			return RETURN_OK;
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
	
	public static void main(String[] args)throws Exception
		{
		/*
		args=new String[]{
				"/commun/data/projects/20120828.AC0KTCACXX.WHOLEGENOME1/align/CD05121/CD05121_recal.bam",
				"/commun/data/users/lindenb/jeter.bam"
				};*/
		new Biostar76892().instanceMain(args);
		}
	}
