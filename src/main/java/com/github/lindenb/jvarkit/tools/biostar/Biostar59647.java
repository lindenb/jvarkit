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
import java.util.Collections;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;

public class Biostar59647 extends AbstractBiostar59647
	{

	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog( Biostar59647.class);

	
	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	static private class MyCommand extends AbstractBiostar59647.AbstractBiostar59647Command
			{
			@Override
			protected Collection<Throwable> call(String inputName) throws Exception
				{
				
			if(super.refFile==null)
				{
				return wrapException("Undefined REF file");
				}
			
		
			IndexedFastaSequenceFile indexedFastaSequenceFile=null;
			SamReader samFileReader=null;
			
			try
				{
				GenomicSequence genomicSequence=null;
				indexedFastaSequenceFile=new IndexedFastaSequenceFile(refFile);
				samFileReader=super.openSamReader(inputName);
				
				
				if(SequenceUtil.areSequenceDictionariesEqual(
						indexedFastaSequenceFile.getSequenceDictionary(),
						samFileReader.getFileHeader().getSequenceDictionary())
						)
					{
					LOG.warn("Not the same sequence dictionaries");
					}
				
				XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
				XMLStreamWriter w= xmlfactory.createXMLStreamWriter(System.out,"UTF-8");
				w.writeStartDocument("UTF-8","1.0");
				w.writeStartElement("sam");
				w.writeAttribute("bam",(inputName==null?"stdin": inputName));
				w.writeAttribute("ref",refFile.getPath());
	
				w.writeComment(getProgramCommandLine());
	
				
				
				
				
	
					
					SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(samFileReader.getFileHeader().getSequenceDictionary());
					SAMRecordIterator iter=samFileReader.iterator();
					while(iter.hasNext())
						{
						SAMRecord rec=iter.next();
						progess.watch(rec);
						final byte readbases[]=rec.getReadBases();
						w.writeStartElement("read");
						
						w.writeStartElement("name");
						w.writeCharacters(rec.getReadName());
						w.writeEndElement();
						
						w.writeStartElement("sequence");
						w.writeCharacters(new String(readbases));
						w.writeEndElement();
						
						w.writeStartElement("flags");
						for(SAMFlag f: SAMFlag.values())
							{
							w.writeAttribute(f.name(),String.valueOf(f.isSet(rec.getFlags())));
							}
											
						w.writeCharacters(String.valueOf(rec.getFlags()));
						w.writeEndElement();//flags
						
						if(!rec.getReadUnmappedFlag())
							{
							w.writeStartElement("qual");
							w.writeCharacters(String.valueOf(rec.getMappingQuality()));
							w.writeEndElement();
							
							w.writeStartElement("chrom");
							w.writeAttribute("index",String.valueOf(rec.getReferenceIndex()));
							w.writeCharacters(rec.getReferenceName());
							w.writeEndElement();
							
							w.writeStartElement("pos");
							w.writeCharacters(String.valueOf(rec.getAlignmentStart()));
							w.writeEndElement();
							
							w.writeStartElement("cigar");
							w.writeCharacters(rec.getCigarString());
							w.writeEndElement();
							}
						
						if(rec.getReadPairedFlag() && !rec.getMateUnmappedFlag())
							{
							w.writeStartElement("mate-chrom");
							w.writeAttribute("index",String.valueOf(rec.getMateReferenceIndex()));
							w.writeCharacters(rec.getMateReferenceName());
							w.writeEndElement();
							
							w.writeStartElement("mate-pos");
							w.writeCharacters(String.valueOf(rec.getMateAlignmentStart()));
							w.writeEndElement();
							}
						
						
						if(!rec.getReadUnmappedFlag() && rec.getCigar()!=null)
							{
							if(genomicSequence==null ||
								!genomicSequence.getChrom().equals(rec.getReferenceName()))
								{
								genomicSequence=new GenomicSequence(indexedFastaSequenceFile, rec.getReferenceName());
								}
							
		
							w.writeStartElement("align");
		
							 int readIndex = 0;
							 int refIndex = rec.getAlignmentStart();
							 				 
							 for (final CigarElement e : rec.getCigar().getCigarElements())
								 {
								 switch (e.getOperator())
									 {
									 case H : break; // ignore hard clips
									 case P : break; // ignore pads
									 case I : //cont.
									 case S :
									 		{
									 		final int length = e.getLength();
									 		for(int i=0;i<length;++i)
									 			{
									 			w.writeEmptyElement(e.getOperator().name());
									 			w.writeAttribute("read-index",String.valueOf(readIndex+1));
									 			if(readIndex>=0 && readIndex< readbases.length)
									 				{
									 				w.writeAttribute("read-base",String.valueOf((char)(readbases[readIndex])));
									 				}
									 			readIndex++;
									 			}
									 		break;
									 		}
									 case N :  //cont. -- reference skip
									 case D :
									 		{
									 		final int length = e.getLength();
									 		for(int i=0;i<length;++i)
									 			{
									 			w.writeEmptyElement(e.getOperator().name());
									 			w.writeAttribute("ref-index",String.valueOf(refIndex));
									 			if(refIndex>=1 && refIndex<= genomicSequence.length())
										 			{
										 			w.writeAttribute("ref-base",String.valueOf(genomicSequence.charAt(refIndex-1)));
										 			}	
									 			refIndex++;
									 			}
									 		break;
									 		}
									 case M :
									 case EQ :
									 case X :
								 			{
									 		final int length = e.getLength();
									 		for(int i=0;i<length;++i)
									 			{
									 			w.writeEmptyElement(e.getOperator().name());
									 			char baseRead='\0';
									 			if(readIndex>=0 && readIndex< readbases.length)
									 				{
										 			baseRead=(char)(rec.getReadBases()[readIndex]);
										 			w.writeAttribute("read-index",String.valueOf(readIndex+1));
										 			w.writeAttribute("read-base",String.valueOf(baseRead));
									 				}
									 			w.writeAttribute("ref-index",String.valueOf(refIndex));
									 			if(refIndex>=1 && refIndex<= genomicSequence.length())
										 			{
									 				char baseRef=genomicSequence.charAt(refIndex-1);
										 			w.writeAttribute("ref-base",String.valueOf(baseRef));
										 			if(Character.toUpperCase(baseRef)!=Character.toUpperCase(baseRead))
									 					{
									 					w.writeAttribute("mismatch","true");
									 					}
										 			}
									 			
									 			
									 			refIndex++;
									 			readIndex++;
									 			}
									 		break;
								 			}
										
									 default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
									 }
								 }
							 w.writeEndElement();
							}
						
						
					
					
					w.writeEndElement();
					}
				iter.close();
				w.writeEndElement();
				w.writeEndDocument();
				w.flush();
				CloserUtil.close(w);
				return Collections.emptyList();
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(samFileReader);
				CloserUtil.close(indexedFastaSequenceFile);
				}
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Biostar59647().instanceMainWithExit(args);

	}

}
