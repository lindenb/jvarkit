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

import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;


public class Biostar84452 extends AbstractBiostar84452
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Biostar84452.class);

	@Override
	public Command createCommand() {
		return new MyCommand();
		}
	
	static private class MyCommand extends AbstractBiostar84452.AbstractBiostar84452Command
		{
		@Override
		protected Collection<Throwable> call(String inputName) throws Exception
			{
			if(super.tag==null || super.tag.length()!=2 || !tag.startsWith("X"))
				{
				return wrapException("Bad tag: expect length=2 && start with 'X'");
				}
				
			
		SAMFileWriter sfw=null;
		SamReader sfr=null;
		try
			{
			sfr = super.openSamReader(inputName);
			SAMFileHeader header=sfr.getFileHeader();
			sfw = openSAMFileWriter(header, true);
			long nChanged=0L;
			SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header);
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec=progress.watch(iter.next());
				if(rec.getReadUnmappedFlag())
					{
					sfw.addAlignment(rec);
					continue;
					}
				
				Cigar cigar=rec.getCigar();
				if(cigar==null)
					{
					sfw.addAlignment(rec);
					continue;
					}
				byte bases[]= rec.getReadBases();
				if(bases==null)
					{
					sfw.addAlignment(rec);
					continue;
					}
				
				ArrayList<CigarElement> L=new ArrayList<CigarElement>();
				ByteArrayOutputStream nseq=new ByteArrayOutputStream();
				ByteArrayOutputStream nqual=new ByteArrayOutputStream();
				
				byte quals[]= rec.getBaseQualities();
				int indexBases=0;
				for(CigarElement ce:cigar.getCigarElements())
					{
					switch(ce.getOperator())
						{
						case S:
							{
							indexBases+=ce.getLength();
							L.add(new CigarElement(ce.getLength(), CigarOperator.H));
							break;
							}
						case H://cont
						case P: //cont
						case N: //cont
						case D:
							{
							L.add(ce);
							break;
							}
							
						case I:
						case EQ:
						case X:
						case M:
							{
							L.add(ce);
							nseq.write(bases,indexBases,ce.getLength());
							if(quals.length!=0) nqual.write(quals,indexBases,ce.getLength());
							indexBases+=ce.getLength(); 
							break;
							}
						default:
							{
							throw new RuntimeException("Unsupported Cigar opertator:"+ce.getOperator());
							}
						}
					
					}
				int ci=0;
				while(ci+1< L.size())
					{
					if(L.get(ci).getOperator().equals(L.get(ci+1).getOperator()))
						{
						L.set(ci, new CigarElement(
								L.get(ci).getLength()+L.get(ci+1).getLength(),
								L.get(ci).getOperator())
								);
						L.remove(ci+1);
						}
					else
						{
						ci++;
						}
					}
				if(indexBases!=bases.length)
					{
					throw new RuntimeException("ERRROR "+rec.getCigarString());
					}
				if(L.size()==cigar.numCigarElements())
					{
					sfw.addAlignment(rec);
					continue;
					}
				++nChanged;
				rec.setAttribute(tag, 1);
				rec.setCigar(new Cigar(L));
				rec.setReadBases(nseq.toByteArray());
				if(quals.length!=0)  rec.setBaseQualities(nqual.toByteArray());
				sfw.addAlignment(rec);
				}
			LOG.info("Num records changed:"+nChanged);
			progress.finish();
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
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new Biostar84452().instanceMainWithExit(args);
		}

	}
