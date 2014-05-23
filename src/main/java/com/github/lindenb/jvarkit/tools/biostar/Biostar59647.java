package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.PrintStream;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;
import com.github.lindenb.jvarkit.util.picard.SamFlag;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;

public class Biostar59647 extends AbstractCommandLineProgram
	{

	private  Biostar59647() {
		}
	
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/Biostar59647";
		}
	
	@Override
	public String getProgramDescription() {
		return "SAM/BAM to XML. See http://www.biostars.org/p/59647/";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -r (reference) Reference file indexed with picard. REQUIRED.");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		File refFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt getopt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=getopt.getopt(args,getGetOptDefault()+ "r:"))!=-1)
			{
			switch(c)
				{
				case 'r': refFile=new File(getopt.getOptArg());break;
				default: 
					{
					switch(handleOtherOptions(c, getopt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default: break;
						}
					}
				}
			}
		
		if(refFile==null)
			{
			error("Undefined REF file");
			return -1;
			}
		File bamFile=null;
		if(getopt.getOptInd()+1!=args.length)
			{
			info("reading from stdin.");
			}
		else
			{
			bamFile=new File(args[getopt.getOptInd()]);
			}
	
		IndexedFastaSequenceFile indexedFastaSequenceFile=null;
		SamReader samFileReader=null;
		
		try
			{
			GenomicSequence genomicSequence=null;
			indexedFastaSequenceFile=new IndexedFastaSequenceFile(refFile);
			samFileReader=null;
			if(bamFile==null)
				{
				samFileReader= SamFileReaderFactory.mewInstance().openStdin();
				}
			else
				{
				samFileReader=SamFileReaderFactory.mewInstance().open(bamFile);
				}
			
			if(SequenceUtil.areSequenceDictionariesEqual(
					indexedFastaSequenceFile.getSequenceDictionary(),
					samFileReader.getFileHeader().getSequenceDictionary())
					)
				{
				warning("Not the same sequence dictionaries");
				}
			
			XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
			XMLStreamWriter w= xmlfactory.createXMLStreamWriter(System.out,"UTF-8");
			w.writeStartDocument("UTF-8","1.0");
			w.writeStartElement("sam");
			w.writeAttribute("ref",(bamFile==null?"stdin": bamFile.getPath()));
			w.writeAttribute("bam", args[1]);

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
					for(SamFlag f: SamFlag.values())
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
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			
			CloserUtil.close(samFileReader);
			CloserUtil.close(indexedFastaSequenceFile);
			}
		return 0;
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Biostar59647().instanceMainWithExit(args);

	}

}
