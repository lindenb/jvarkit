package com.github.lindenb.jvarkit.tools.blast2sam;

import gov.nih.nlm.ncbi.blast.Hit;
import gov.nih.nlm.ncbi.blast.Hsp;
import gov.nih.nlm.ncbi.blast.Iteration;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.bio.blast.BlastHspAlignment;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;

public class BlastToSam extends AbstractCommandLineProgram
	{
	private SAMSequenceDictionary dictionary;
	private Unmarshaller unmarshaller;
	private int EXPECTED_SIZE=500;
	//fool javac
	@SuppressWarnings("unused")
	private final static gov.nih.nlm.ncbi.blast.ObjectFactory _foolJavac=null;
	
	private static class SequenceIteration
		{
		//Iteration iteration;
		//String queryDef;
		List<SAMRecord> records=new ArrayList<SAMRecord>();
		}
	
	private class Paired
		implements Comparable<Paired>
		{
		SAMRecord rec1;
		SAMRecord rec2;
		
		void completeFlags()
			{
			rec1.setFirstOfPairFlag(true);
			rec1.setSecondOfPairFlag(false);
			
			rec2.setFirstOfPairFlag(false);
			rec2.setSecondOfPairFlag(true);
			

			rec1.setReadPairedFlag(true);
			rec2.setReadPairedFlag(true);
			
			rec1.setProperPairFlag(false);
			rec2.setProperPairFlag(false);
			
			for(int i=0;i< 2;++i)
				{
				SAMRecord recA=(i==0?rec1:rec2);
				if(recA.getReadUnmappedFlag())
					{
					recA.setReferenceIndex(-1);
					}
				}
			
			for(int i=0;i< 2;++i)
				{
				SAMRecord recA=(i==0?rec1:rec2);
				SAMRecord recB=(i==0?rec2:rec1);
				recA.setMateUnmappedFlag(recB.getReadUnmappedFlag());
				
				
				if(!recB.getReadUnmappedFlag())
					{
					recA.setMateReferenceName(recB.getReferenceName());
					recA.setMateReferenceIndex(recB.getReferenceIndex());
					recA.setMateNegativeStrandFlag(recB.getReadNegativeStrandFlag());
					recA.setMateAlignmentStart(recB.getAlignmentStart());
					}
				else
					{
					recA.setMateReferenceIndex(-1);
					
					}
				}	
			if( !rec1.getReadUnmappedFlag() &&
				!rec2.getReadUnmappedFlag() && 
				rec1.getReadNegativeStrandFlag()!=rec2.getReadNegativeStrandFlag() &&
				rec1.getReferenceName().equals(rec2.getReferenceName())		 
				)
				{
				int len=rec1.getAlignmentStart()-rec2.getAlignmentEnd();
				if(Math.abs(len) <= EXPECTED_SIZE)
					{
					rec1.setProperPairFlag(true);
					rec2.setProperPairFlag(true);
					}
				rec1.setInferredInsertSize(-len);
				rec2.setInferredInsertSize(len);
				
				}
			}
		
		int score()
			{
			int v=0;
			v+=rec1.getReadUnmappedFlag()?0:1;
			v+=rec2.getReadUnmappedFlag()?0:1;
			if( !rec1.getReadUnmappedFlag() && !rec2.getReadUnmappedFlag()  )		 
				{
				if(rec1.getReferenceName().equals(rec2.getReferenceName()))
					{
					v+=2;
					if(rec1.getReadNegativeStrandFlag()!=rec2.getReadNegativeStrandFlag())
						{
						v+=4;
						int len=rec1.getAlignmentStart()-rec2.getAlignmentEnd();
						if(Math.abs(len) <= EXPECTED_SIZE)
							{
							v+=8;
							}
						
						}
					}
				}
			return v;
			}
		
		@Override
		public int compareTo(Paired p)
			{
			return p.score()-score();//greater is best
			}
		}
	
	private BlastToSam()
		{
		
		}
	
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/Blast2Sam";
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
	
	private void dumpSingle(SAMFileWriter w,SequenceIteration si)
		{
		boolean first=true;
		for(SAMRecord rec:si.records)
			{
			rec.setNotPrimaryAlignmentFlag(!first);
			first=false;
			w.addAlignment(rec);
			}
		si.records.clear();
		}
	
	private void run_single(
			SAMFileWriter w,
			XMLEventReader r,
			SAMFileHeader header
			)
			throws XMLStreamException,JAXBException
		{
		List<Iteration> stack=new ArrayList<Iteration>();
		String prev=null;
		for(;;)
			{
			Iteration iter1=peekIteration(r);
			if(iter1==null || !(iter1.getIterationQueryDef().equals(prev)))
				{
				SequenceIteration si=convertIterationToSequenceIteration(stack,header);
				dumpSingle(w,si);
				if(iter1==null) break;
				stack.clear();
				prev=iter1.getIterationQueryDef();
				}
			stack.add(iter1);
			}
		}
	
	private SequenceIteration convertIterationToSequenceIteration(
			List<Iteration> stack,
			final SAMFileHeader header
			)
			throws XMLStreamException,JAXBException
			{
			SequenceIteration sequenceIteration=new SequenceIteration(); 
			if(stack.isEmpty()) return sequenceIteration;
			
			SAMReadGroupRecord rg1=header.getReadGroup("g1");
			//sequenceIteration.iteration=iter1;
			
			SAMRecordFactory samRecordFactory=new DefaultSAMRecordFactory();

			
			
			final StringBuilder readContent=new StringBuilder();
			final int iterLength=Integer.parseInt(stack.get(0).getIterationQueryLen());
			
			for(Iteration iter1:stack)
				{
				for(Hit hit: iter1.getIterationHits().getHit())
					{
					
					for(Hsp hsp: hit.getHitHsps().getHsp())
						{
						
						for(BlastHspAlignment.Align a:new BlastHspAlignment(hsp))
							{
							char c=a.getQueryChar();
							if(!Character.isLetter(c)) continue;
							int queryIndex0=a.getQueryIndex1()-1;
							while(readContent.length()<=queryIndex0) readContent.append('N');
							if(readContent.charAt(queryIndex0)=='N')
								{
								readContent.setCharAt(queryIndex0, c);
								}
							else if(readContent.charAt(queryIndex0)!=c)
								{
								throw new IllegalStateException();
								}
							}
						}
					}
				}
			
			
			for(Iteration iter1:stack)
				{
				for(Hit hit: iter1.getIterationHits().getHit())
					{
					for(Hsp hsp: hit.getHitHsps().getHsp())
						{
						SAMRecord rec=samRecordFactory.createSAMRecord(header);
						rec.setReadUnmappedFlag(false);
						rec.setReadName(iter1.getIterationQueryDef());
						if( hit.getHitAccession()!=null &&
							!hit.getHitAccession().trim().isEmpty() &&
							this.dictionary.getSequence(hit.getHitAccession())!=null
							)
							{
							rec.setReferenceName(hit.getHitAccession());
							}
						else
							{
							rec.setReferenceName(hit.getHitDef());
							}
						SAMSequenceRecord ssr=this.dictionary.getSequence(hit.getHitDef());
						if(ssr==null)
							{
							warning("Hit is not in SAMDictionary "+hit.getHitDef());
							rec.setReferenceIndex(-1);
							}
						else
							{
							rec.setReferenceIndex(ssr.getSequenceIndex());
							}
						
						BlastHspAlignment blastHspAlignment=new BlastHspAlignment(hsp);
						rec.setReadNegativeStrandFlag(blastHspAlignment.isPlusMinus());
	
						
						List<CigarOperator> cigarL=new ArrayList<CigarOperator>();
						for(BlastHspAlignment.Align a:blastHspAlignment)
							{
							//System.err.println("##"+a);
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
	
						
						if(cigarL.size()!=hsp.getHspMidline().length())
							{
							throw new IllegalStateException("Boumm");
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
						/* soft clip */ 
						if(blastHspAlignment.getQueryTo1()< readContent.length())
							{
							cigarE.add(new CigarElement(
									(readContent.length()-blastHspAlignment.getQueryTo1()),
									CigarOperator.S 
									));
							}
						/* hard clip */
						if(readContent.length() < iterLength)
							{
							cigarE.add(new CigarElement(
									(iterLength-readContent.length()),
									CigarOperator.H
									));
							}
						
						
						rec.setCigar(cigarE);
						rec.setMappingQuality(40);
						rec.setAlignmentStart(Math.min(blastHspAlignment.getHitFrom1(),blastHspAlignment.getHitTo1()));
						rec.setAttribute("BB", Float.parseFloat(hsp.getHspBitScore()));
						rec.setAttribute("BE", Float.parseFloat(hsp.getHspEvalue()));
						rec.setAttribute("BS", Float.parseFloat(hsp.getHspScore()));
						rec.setAttribute("NM", Integer.parseInt(hsp.getHspGaps()));
						rec.setAttribute("RG", rg1.getId());
						// setAlignmentEnd not supported in SAM API
						//rec.setAlignmentEnd(Math.max(blastHspAlignment.getHitFrom1(),blastHspAlignment.getHitTo1())); 
						sequenceIteration.records.add(rec);
						}
					}
				}
			
			if(readContent.length()==0)
				{
				readContent.append('N');
				}
			
			byte readBases[]=readContent.toString().getBytes();
			char readQuals[]=new char[readBases.length];
			
			 
			
			for(int i=0;i< readBases.length;++i)
				{
				readQuals[i]=(readBases[i]=='N'?'#':'J');
				}

			
			
			if(sequenceIteration.records.isEmpty())
				{
				SAMRecord rec=samRecordFactory.createSAMRecord(header);
				rec.setReadName(stack.get(0).getIterationQueryDef());
				rec.setReadUnmappedFlag(true);
				rec.setAttribute("RG", rg1.getId());
				sequenceIteration.records.add(rec);
				}
			
			
				
			
			for(SAMRecord rec:sequenceIteration.records)
				{
				rec.setReadString(new String(readBases));
				rec.setReadBases(readBases);
				rec.setBaseQualityString(new String(readQuals,0,readQuals.length));
				rec.setBaseQualities(htsjdk.samtools.SAMUtils.fastqToPhred(new String(readQuals,0,readQuals.length)));
				}
		return sequenceIteration;
		}
	
	private static SAMRecord cloneSAMRecord(final SAMRecord rec)
		{
		try {
			return (SAMRecord)rec.clone();
			}
		catch (Exception e)
			{
			throw new RuntimeException("Cannot clone a SAMRecord ?",e);
			}
		}
	
	private void dumpPaired(SAMFileWriter w,SequenceIteration si1,SequenceIteration si2)
		{
		if(si1.records.isEmpty()) return;
		
		SequenceIteration siL[]=new SequenceIteration[]{si1,si2};
		for(SequenceIteration si:siL)
			{
			for(SAMRecord rec:si.records)
				{
				rec.setReadPairedFlag(true);
				rec.setMateUnmappedFlag(true);
				}
			}
		List<Paired> paired=new ArrayList<Paired>();
		int x=0;
		while(x < si1.records.size())
			{
			SAMRecord rec1=si1.records.get(x);
			int y=0;
			
			while(y < si2.records.size())
				{
				SAMRecord rec2=si2.records.get(y);
				
				Paired pair=new Paired();
				
				pair.rec1=cloneSAMRecord(rec1);
				pair.rec2=cloneSAMRecord(rec2);
				
				if(pair.rec1.getReadUnmappedFlag() && pair.rec2.getReadUnmappedFlag())
					{
					++y;
					continue;
					}	
				
				if(!paired.isEmpty() )
					{
					int cmp= pair.compareTo(paired.get(0)) ;
					if(cmp<0)
						{
						paired.clear();
						}
					else if(cmp>0)
						{
						++y;
						continue;
						}
					}
				paired.add(pair);
				++y;
				}
			++x;
			}
		
		if(paired.isEmpty())
			{
			Paired pair=new Paired();
			pair.rec1=cloneSAMRecord(si1.records.get(0));
			pair.rec2=cloneSAMRecord(si2.records.get(0));
			paired.add(pair);
			}
		
		for(int i=0;i< paired.size();++i)
			{
			Paired pair=paired.get(i);
			pair.completeFlags();
			if(!pair.rec1.getReadUnmappedFlag())
				{
				pair.rec1.setNotPrimaryAlignmentFlag(i!=0);
				}
			if(!pair.rec2.getReadUnmappedFlag())
				{
				pair.rec2.setNotPrimaryAlignmentFlag(i!=0);
				}
			w.addAlignment(pair.rec1);
			w.addAlignment(pair.rec2);
			}
		si1.records.clear();
		si2.records.clear();
		}
	
	private void run_paired(
			SAMFileWriter w,
			XMLEventReader r,
			SAMFileHeader header
			)
			throws XMLStreamException,JAXBException
		{
		List<Iteration> stack1=new ArrayList<Iteration>();
		Iteration iter=null;
		for(;;)
			{
			String prev_name=null;
			if( iter==null)
				{
				iter=peekIteration(r);
				if(iter==null) break;
				}
			stack1.add(iter);
			List<Iteration> stack2=new ArrayList<Iteration>();
			prev_name=iter.getIterationQueryDef();
			
			//pileup first of pair
			for(;;)
				{
				iter=peekIteration(r);
				if(iter==null)
					{
					throw new RuntimeException("Illegal number of read forward/reverse");
					}
				else if(iter.getIterationQueryDef().equals(prev_name))
					{
					stack1.add(iter);
					}
				else
					{
					stack2.add(iter);
					prev_name=iter.getIterationQueryDef();
					break;
					}
				}
			
			//pileup second of pair
			for(;;)
				{
				iter=peekIteration(r);
				if(iter==null || !iter.getIterationQueryDef().equals(prev_name))
					{
					SequenceIteration si1=convertIterationToSequenceIteration(stack1, header);
					SequenceIteration si2=convertIterationToSequenceIteration(stack2, header);
					dumpPaired(w,si1,si2);
					stack1.clear();
					stack2.clear();
					break;
					}
				else
					{
					stack2.add(iter);
					}
				}
			if(iter==null) break;
			}
		
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -r (file)  fasta sequence file indexed with picard. Required.");
		out.println(" -o (file.bam) filename out . Default: SAM stdout");
		out.println(" -p (int expected size) input is an interleaved list of sequences forward and reverse (paired-ends)");
		super.printOptions(out);
		}

	
	@Override
	public int doWork(String[] args)
		{
		boolean interleaved_input=false;
		int maxRecordsInRam=10000;
		File fileout=null;
		String faidx=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt getopt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=getopt.getopt(args, super.getGetOptDefault()+ "r:o:p:"))!=-1)
			{
			switch(c)
				{	
				case 'p':
					interleaved_input=true;
						EXPECTED_SIZE=Integer.parseInt(getopt.getOptArg());
						break;
				case 'r': faidx=getopt.getOptArg();break;
				case 'o': fileout=new File(getopt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, getopt, null))
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
			
			
			SAMProgramRecord prg2=header.createProgramRecord();
			fillHeader(rx,prg2);
			SAMProgramRecord prg1=header.createProgramRecord();
			prg1.setCommandLine(getProgramCommandLine());
			prg1.setProgramVersion(getVersion());
			prg1.setProgramName(getProgramName());
			prg1.setPreviousProgramGroupId(prg2.getId());
			SAMReadGroupRecord rg1=new SAMReadGroupRecord("g1");
			rg1.setLibrary("blast");
			rg1.setSample("blast");
			rg1.setDescription("blast");
			header.addReadGroup(rg1);
			
			

			if(fileout==null)
				{
				sfw=sfwf.makeSAMWriter(header, false, System.out);
				}
			else
				{
				sfw=sfwf.makeSAMOrBAMWriter(header, false, fileout);
				}

			
			if(interleaved_input)
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
