package com.github.lindenb.jvarkit.tools.blast2sam;

import gov.nih.nlm.ncbi.blast.Hit;
import gov.nih.nlm.ncbi.blast.Hsp;
import gov.nih.nlm.ncbi.blast.Iteration;
import gov.nih.nlm.ncbi.blast.IterationHits;

import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import net.sf.picard.PicardException;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.CigarUtil;
import net.sf.samtools.DefaultSAMRecordFactory;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordFactory;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.bio.blast.BlastHspAlignment;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;

public class BlastToSam extends AbstractCommandLineProgram
	{
	private SAMSequenceDictionary dictionary;
	private Unmarshaller unmarshaller;
	//fool javac
	private final static gov.nih.nlm.ncbi.blast.ObjectFactory _foolJavac=null;
	
	private BlastToSam()
		{
		
		}
	
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/BlastToSam";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -r (file)  fasta sequence file indexed with picard. Required.");
		super.printOptions(out);
		}

	private Iteration peekIteration(XMLEventReader r) throws XMLStreamException,JAXBException
		{
		while(r.hasNext())
			{
			XMLEvent evt=r.peek();
		
			if(!(evt.isStartElement() && evt.asStartElement().getName().getLocalPart().equals("Iteration")))
				{
				r.next();
				continue;
				}
			return this.unmarshaller.unmarshal(r, Iteration.class).getValue();
			}
		return null;
		}
	
	private void fillHeader(XMLEventReader r,SAMProgramRecord prog) throws XMLStreamException,JAXBException
		{
		while(r.hasNext())
			{
			XMLEvent evt=r.peek();
		
			if(!(evt.isStartElement()))
				{
				r.next();
				continue;
				}
			StartElement E=evt.asStartElement();
			String name=E.getName().getLocalPart();
			if(name.equals("BlastOutput_iterations")) break;
			r.next();
			if(name.equals("BlastOutput_program"))
				{
				prog.setProgramName(r.getElementText());
				}
			else if(name.equals("BlastOutput_version"))
				{
				prog.setProgramVersion(r.getElementText().replace(' ', '_'));
				}
			}
		}
	
	private void run_single(
			SAMFileWriter w,
			XMLEventReader r,
			SAMFileHeader header
			)
			throws XMLStreamException,JAXBException
		{
		
		for(;;)
			{
			SAMRecordFactory samRecordFactory=new DefaultSAMRecordFactory();

			Iteration iter1=peekIteration(r);
			if(iter1==null) return;
			if(iter1.getIterationHits().getHit().isEmpty())
				{
				SAMRecord rec=samRecordFactory.createSAMRecord(header);
				rec.setReadName(iter1.getIterationQueryDef());
				rec.setReadUnmappedFlag(true);
				w.addAlignment(rec);
				}
			else
				{
				boolean first=true;
				
				StringBuilder readContent=new StringBuilder();
				int iterLength=Integer.parseInt(iter1.getIterationQueryLen());
				while(readContent.length()<iterLength ) readContent.append('N');
				
				for(Hit hit: iter1.getIterationHits().getHit())
					{
					
					for(Hsp hsp: hit.getHitHsps().getHsp())
						{
						
						for(BlastHspAlignment.Align a:new BlastHspAlignment(hsp))
							{
							char c=a.getQueryChar();
							if(!Character.isLetter(c)) continue;
							int read1=a.getQueryIndex1();
							while(readContent.length()<read1) readContent.append('N');
							if(readContent.charAt(read1-1)=='N')
								{
								readContent.setCharAt(read1-1, c);
								}
							else if(readContent.charAt(read1-1)!=c)
								{
								throw new IllegalStateException();
								}
							}
						}
					}
				
				byte readBases[]=readContent.toString().getBytes();
				char readQuals[]=new char[readBases.length];
				
				 
				
				for(int i=0;i< readBases.length;++i)
					{
					readQuals[i]=(readBases[i]=='N'?'#':'J');
					}
				
				for(Hit hit: iter1.getIterationHits().getHit())
					{
					for(Hsp hsp: hit.getHitHsps().getHsp())
						{
						SAMRecord rec=samRecordFactory.createSAMRecord(header);
						rec.setReadName(iter1.getIterationQueryDef());
						rec.setNotPrimaryAlignmentFlag(!first);
						if(first) first=false;
						if(hit.getHitAccession()!=null && !hit.getHitAccession().trim().isEmpty())
							{
							rec.setReferenceName(hit.getHitAccession());
							}
						else
							{
							rec.setReferenceName(hit.getHitDef());
							}
						rec.setReadString(new String(readBases));
						rec.setReadBases(readBases);
						rec.setBaseQualityString(new String(readQuals,0,readQuals.length));
						rec.setBaseQualities(net.sf.samtools.SAMUtils.fastqToPhred(new String(readQuals,0,readQuals.length)));

						
						BlastHspAlignment blastHspAlignment=new BlastHspAlignment(hsp);
						List<CigarOperator> cigarL=new ArrayList<CigarOperator>();
						for(BlastHspAlignment.Align a:blastHspAlignment)
							{
							if(a.getMidChar()=='|')
								{
								cigarL.add(CigarOperator.EQ);
								}
							else if(a.getMidChar()==':')
								{
								cigarL.add(CigarOperator.M);
								}
							else if(a.getHitChar()=='-')
								{
								cigarL.add(CigarOperator.I);
								}
							else if(a.getQueryChar()=='-')
								{
								cigarL.add(CigarOperator.D);
								}
							else
								{
								cigarL.add(CigarOperator.X);
								}
							}
						
						
						
						Cigar cigarE=new Cigar();
						
						if(blastHspAlignment.getQueryFrom1()>1)
							{
							cigarE.add(new CigarElement(
									blastHspAlignment.getQueryFrom1()-1,
									CigarOperator.S
									));
							}
						int x=0;
						while(x< cigarL.size())
							{
							int y=x+1;
							while(y< cigarL.size() && cigarL.get(x)==cigarL.get(y))
								{
								++y;
								}
							cigarE.add(new CigarElement(y-x, cigarL.get(x)));
							x=y;
							}
						if(blastHspAlignment.getQueryTo1()< iterLength)
							{
							cigarE.add(new CigarElement(
									iterLength-blastHspAlignment.getQueryTo1(),
									CigarOperator.S
									));
							}
						
						
						rec.setCigar(cigarE);
						rec.setMappingQuality(40);
						rec.setAlignmentStart(Math.min(blastHspAlignment.getHitFrom1(),blastHspAlignment.getHitTo1()));
						// setAlignmentEnd not supported in SAM API
						//rec.setAlignmentEnd(Math.max(blastHspAlignment.getHitFrom1(),blastHspAlignment.getHitTo1())); 
						w.addAlignment(rec);
						}
					}
				}
			}
		}
	
	private void run_paired(
			SAMFileWriter w,
			XMLEventReader r,
			SAMFileHeader header
			)
			throws XMLStreamException,JAXBException
		{
		for(;;)
			{
			Iteration read1=peekIteration(r);
			if(read1==null) return;
			Iteration read2=peekIteration(r);
			if(read2==null) throw new PicardException("Illegal number of read forward/reverse");
			
			}
		}
	
	
	
	@Override
	public int doWork(String[] args)
		{
		boolean interlaced_input=false;
		int maxRecordsInRam=10000;
		File fileout=null;
		String faidx=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt getopt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=getopt.getopt(args, super.getGetOptDefault()+ "r:o:"))!=-1)
			{
			switch(c)
				{
				case 'r': faidx=getopt.getOptArg();break;
				case 'o': fileout=new File(getopt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, getopt))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		if(faidx==null)
			{
			error("Indexed fasta file missing.");
			return -1;
			}
		SAMFileWriter sfw=null;
		XMLEventReader rx=null;
		SAMFileWriterFactory sfwf=new SAMFileWriterFactory();
		sfwf.setCreateIndex(false);
		sfwf.setMaxRecordsInRam(maxRecordsInRam);
		sfwf.setCreateMd5File(false);
		sfwf.setUseAsyncIo(false);
		SAMFileHeader header=new SAMFileHeader();
		try
			{
			info("opening "+faidx);
			this.dictionary=new SAMSequenceDictionaryFactory().load(new File(faidx));
			header.setSortOrder(SortOrder.unsorted);
			header.setSequenceDictionary(this.dictionary);
			
			
			if(fileout==null)
				{
				sfw=sfwf.makeSAMWriter(header, false, System.out);
				}
			else
				{
				sfw=sfwf.makeSAMOrBAMWriter(header, false, fileout);
				}
			JAXBContext jc = JAXBContext.newInstance("gov.nih.nlm.ncbi.blast");
			this.unmarshaller=jc.createUnmarshaller();
			XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
			xmlInputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.FALSE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, Boolean.FALSE);
			xmlInputFactory.setXMLResolver(new XMLResolver()
				{
				@Override
				public Object resolveEntity(String arg0, String arg1, String arg2,
						String arg3) throws XMLStreamException
					{
					info("resolveEntity:" +arg0+"/"+arg1+"/"+arg2);
					return null;
					}
				});
			if(getopt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				rx=xmlInputFactory.createXMLEventReader(System.in);
				}
			else if(getopt.getOptInd()+1==args.length)
				{
				info("Reading from "+args[getopt.getOptInd()]);
				rx=xmlInputFactory.createXMLEventReader(IOUtils.openURIForBufferedReading(args[getopt.getOptInd()]));
				}
			else
				{
				error("Illegal number of args");
				return -1;
				}
			
			SAMProgramRecord prg=header.createProgramRecord();
			fillHeader(rx,prg);
			
			if(interlaced_input)
				{
				run_paired(sfw,rx,header);
				}
			else
				{
				run_single(sfw,rx,header);
				}
			
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}	
		finally
			{
			CloserUtil.close(sfw);
			CloserUtil.close(rx);
			}
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new BlastToSam().instanceMainWithExit(args);
		}

	}
