package com.github.lindenb.jvarkit.tools.cmpbams;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.Marshaller;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.tools.cmpbams.entities.BamRecord;
import com.github.lindenb.jvarkit.tools.cmpbams.entities.Records;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.IntervalUtils;

import com.github.lindenb.jvarkit.util.picard.cmdline.Option;
import com.github.lindenb.jvarkit.util.picard.cmdline.StandardOptionDefinitions;
import com.github.lindenb.jvarkit.util.picard.cmdline.Usage;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;


public class CompareBams3  extends AbstractCommandLineProgram
	{
	private static final Logger LOG=Logger.getLogger(CompareBams3.class.getName());
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"Compare two or more BAM files, generate a XML report.";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="BAM files to process.",
    		minElements=2,optional=false)
	public List<File> IN=new ArrayList<File>();
    @Option(shortName= "L", doc="restrict to that region (chr:start-end)",optional=true)
	public String REGION=null;
	
    private class ReadComparator
	implements Comparator<Match>
		{
		@Override
		public int compare(Match m0, Match m1)
			{
			int i=m0.readName.compareTo(m1.readName);
			if(i!=0) return i;
			i=m0.num_in_pair-m1.num_in_pair;
			return i;
			}
		}
	
	
	private class MatchCodec
		extends AbstractDataCodec<Match>
		{
		@Override
		public MatchCodec clone()
			{
			return new MatchCodec();
			}
		@Override
		public Match decode(DataInputStream dis) throws IOException
			{
			Match m=new Match();
			try {
				m.readName=dis.readUTF();
				}
			catch(IOException err)
				{
				return null;
				}
			m.bamIndex=dis.readInt();
			m.tid=dis.readInt();
			m.pos=dis.readInt();
			m.num_in_pair=dis.readByte();
			m.flag=dis.readInt();
			m.cigar=dis.readUTF();
			
			m.rnext=dis.readInt();
			m.pnext=dis.readInt();
			m.tlen=dis.readInt();
			m.mapq=dis.readInt();
			
			return m;
			}
		@Override
		public void encode(DataOutputStream dos, Match match)
				throws IOException
			{
			dos.writeUTF(match.readName);
			dos.writeInt(match.bamIndex);
			dos.writeInt(match.tid);
			dos.writeInt(match.pos);
			dos.writeByte(match.num_in_pair);
			dos.writeInt(match.flag);
			dos.writeUTF(match.cigar);
			dos.writeInt(match.rnext);
			dos.writeInt(match.pnext);
			dos.writeInt(match.tlen);
			dos.writeInt(match.mapq);
			
			}
		
		}
	
	private static class Match
		{
		String readName="";
		byte num_in_pair=0;
		int tid=-1;
		int bamIndex=-1;
		int pos=-1;
		int flag=0;
		String cigar="";
		int rnext=-1;
		int pnext=-1;
		int tlen=0;
		int mapq=0;
		}
	
	
	
	@Override
	protected int doWork()
		{
		SAMFileReader samFileReader=null;
		try
			{
			if(this.IN.size() <2)
				{
				System.err.println("Need more bams please");
				return -1;
				}
			
			
			final ReadComparator matchComparator=new ReadComparator();
			SortingCollection<Match> database=SortingCollection.newInstance(
					Match.class,
					new MatchCodec(),
					matchComparator,
					super.MAX_RECORDS_IN_RAM
					);
			database.setDestructiveIteration(true);
	
			List<SAMSequenceDictionary> sequenceDictionaries=new ArrayList<SAMSequenceDictionary>(this.IN.size());
			
			for(int currentSamFileIndex=0;
					currentSamFileIndex<this.IN.size();
					currentSamFileIndex++ )
				{
				long nReads=0L;
				File samFile=this.IN.get(currentSamFileIndex);
				LOG.info("Opening "+samFile);
				samFileReader=new SAMFileReader(samFile);
				samFileReader.setValidationStringency(super.VALIDATION_STRINGENCY);
				SAMSequenceDictionary dict=samFileReader.getFileHeader().getSequenceDictionary();
				sequenceDictionaries.add(dict);
				
				Interval interval=null;
				if(REGION!=null)
					{
					interval=IntervalUtils.parseOne(dict, REGION);
					if(interval==null)
						{
						System.err.println("Cannot parse "+REGION+" (bad syntax or not in dictionary");
						return -1;
						}
					}
				
				Iterator<SAMRecord> iter=null;
				if(interval==null)
					{
					iter=samFileReader.iterator();
					}
				else
					{
					iter=samFileReader.queryOverlapping(interval.getSequence(), interval.getStart(), interval.getEnd());
					}
				
				while(iter.hasNext() )
					{
					if(nReads++%10000000==0) LOG.info("in "+samFile+" count:"+nReads);
					
					
					
					SAMRecord rec=iter.next();
					
					Match m=new Match();
					if(rec.getReadPairedFlag())
						{
						m.num_in_pair=(byte)(rec.getFirstOfPairFlag()?1:2);
						
						if(!rec.getMateUnmappedFlag())
							{
							m.rnext=rec.getMateReferenceIndex();
							m.pnext=rec.getMateAlignmentStart();
							}
						else
							{
							m.rnext=-1;
							m.pnext=-1;
							}
						
						}
					else
						{
						m.num_in_pair=0;
						}
					m.readName=rec.getReadName();
					m.bamIndex=currentSamFileIndex;
					m.flag=rec.getFlags();
					if(rec.getReadUnmappedFlag())
						{
						m.tid=-1;
						m.pos=-1;
						m.mapq=0;
						m.tlen=0;
						m.cigar="";
						}
					else
						{
						m.tid=rec.getReferenceIndex();
						m.pos=rec.getAlignmentStart();
						m.mapq=rec.getMappingQuality();
						m.cigar=rec.getCigarString();
						if(m.cigar==null) m.cigar="";
						if(m.rnext!=-1)//mate mapped
							{
							m.tlen=rec.getInferredInsertSize();
							}
						}
					
					database.add(m);
					}
				samFileReader.close();
				samFileReader=null;
				LOG.info("Close "+samFile);
				}
			database.doneAdding();
			LOG.info("Writing results....");
			
			XMLOutputFactory xmlOutputFactory=XMLOutputFactory.newFactory();
			XMLStreamWriter w=xmlOutputFactory.createXMLStreamWriter(System.out, "UTF-8");
			w.writeStartDocument("UTF-8", "1.0");
			w.writeStartElement("comparebams");
			
			w.writeStartElement("header");
			w.writeStartElement("files");
			for(int x=0;x<this.IN.size();++x)
				{
				w.writeEmptyElement("file");
				w.writeAttribute("path", this.IN.get(x).getPath());
				w.writeAttribute("index", String.valueOf(x+1));
				}
			w.writeEndElement();//files
			w.writeEndElement();//header
			
			w.writeStartElement("body");
			w.writeCharacters("\n");
			
			
			
			/* create an array of set<Match> */
			List<Match> row=new ArrayList<CompareBams3.Match>(this.IN.size());
			
			JAXBContext jc = JAXBContext.newInstance(BamRecord.class,Records.class);
			Marshaller marshaller=jc.createMarshaller();
			marshaller.setProperty(Marshaller.JAXB_FRAGMENT, Boolean.TRUE);
			CloseableIterator<Match> iter=database.iterator();
			for(;;)
				{
				Match nextMatch = null;
				if(iter.hasNext())
					{
					nextMatch = iter.next();
					}
				if(nextMatch==null && row.isEmpty()) break;
				if(!row.isEmpty() &&  (nextMatch==null || matchComparator.compare(nextMatch, row.get(0))!=0))
					{
					Match first=row.get(0);
					Records records=new Records();
					records.setName(first.readName);
					records.setSide(first.num_in_pair);
					for(Match m:row)
						{
						BamRecord br=new BamRecord();
						br.setFileIndex(m.bamIndex+1);
						if(m.tid!=-1)
							{
							br.setChrom(sequenceDictionaries.get(m.bamIndex).getSequence(m.tid).getSequenceName());
							br.setPos(m.pos);
							if(m.rnext!=-1)
								{
								br.setTlen(m.tlen);
								}
							br.setMapq(m.mapq);
							if(!m.cigar.isEmpty()) br.setCigar(m.cigar);
							}
						if(m.rnext!=-1)
							{
							br.setRnext(sequenceDictionaries.get(m.bamIndex).getSequence(m.rnext).getSequenceName());
							br.setPnext(m.pnext);
							}
						br.setFlag(m.flag);
						
						
						records.getRecord().add(br);
						
						}
					
					marshaller.marshal(new JAXBElement<Records>(
							new QName("records"),
							Records.class, records), w);
					w.writeCharacters("\n");
					if(nextMatch==null) break;
					row.clear();
					}
				row.add(nextMatch);
				}
			
			iter.close();
			
			w.writeEndElement();//body
			w.writeEndElement();//cmpbames
			w.writeEndDocument();
			w.flush();
			w.close();
			
			}
		catch(Exception err)
			{
			err.printStackTrace();
			return -1;
			}
		finally
			{
			if(samFileReader!=null) samFileReader.close();
			}
		return 0;
		}
		
	public static void main(String[] args) throws Exception
		{
		new CompareBams3().instanceMainWithExit(args);
		}
}
