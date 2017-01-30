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

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;


import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;

public class Biostar234081 extends AbstractBiostar234081
	{
	
@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
	SamReader in =null;
	SAMFileWriter w=null;  
	SAMRecordIterator iter=null;  
	try {
		in = super.openSamReader(inputName);
		w = super.openSAMFileWriter(in.getFileHeader(),true);
		iter=in.iterator();
		while(iter.hasNext())
			{
			SAMRecord rec = iter.next();
			if(!rec.getReadUnmappedFlag() && rec.getCigar()!=null)
				{
				final Cigar cigar = rec.getCigar();
				final List<CigarElement> elements = new ArrayList<>(cigar.numCigarElements());

				for(int i=0;i< cigar.numCigarElements();++i)
					{
					CigarElement ce =cigar.getCigarElement(i);
					switch(ce.getOperator())
						{
						case X:
						case EQ: ce =new CigarElement(ce.getLength(), CigarOperator.M);break;
						default:break;
						}
					if(!elements.isEmpty() && elements.get(elements.size()-1).getOperator()==ce.getOperator())
						{
						elements.set(elements.size()-1,
								new CigarElement(
										ce.getLength()+elements.get(elements.size()-1).getLength(),
										ce.getOperator())
								);
						}
					else
						{
						elements.add(ce);
						}
					}
				rec.setCigar(new Cigar(elements));
				}
			w.addAlignment(rec);
			}
		return RETURN_OK;
	} catch (final Exception err) {
		return wrapException(err);	
	}	
	finally {
		CloserUtil.close(iter);
		CloserUtil.close(in);
		CloserUtil.close(w);
	
	}
	
	
		}	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new Biostar234081().instanceMainWithExit(args);
		}

	}
