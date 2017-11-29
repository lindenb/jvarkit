/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.function.IntFunction;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.misc.ConcatSam;
import com.github.lindenb.jvarkit.util.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamFilterParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
/**

BEGIN_DOC

END_DOC
*/
@Program(name="samtranslocations",
	description="Explore translocations between two chromosomes",
	keywords={"sam","bam"},
	generate_doc=false
	)
public class SamTranslocations extends Launcher {
	private static final Logger LOG = Logger.build(SamTranslocations.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"--region","--interval"},description="Limit analysis to this interval. "+IntervalParser.OPT_DESC)
	private String region_str=null;
	@Parameter(names={"--filter"},description=SamFilterParser.FILTER_DESCRIPTION,converter=SamFilterParser.StringConverter.class)
	private SamRecordFilter samRecordFilter = SamFilterParser.buildDefault();
	@Parameter(names={"--groupby"},description="Group Reads by. "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition samRecordPartition = SAMRecordPartition.sample;
	@Parameter(names={"-m","--min-events"},description="Minimal number of events for printing a result")
	private int min_number_of_event = 1;
	@Parameter(names={"-r","--round"},description="Round locations to LOC=LOC-LOC%round")
	private int round_loc = 150;
	@Parameter(names={"-x","--xml"},description="XML Output")
	private boolean xml_output=false;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();

	private static class Event
		{
		String partition;
		int tid1;
		int pos1;
	    int tid2;
	    int pos2;
	    long id;
		Event() {
			
			}
		
		int compare2(final Event o) {
			int i= partition.compareTo(o.partition);
			if(i!=0) return i;
			i = tid1 - o.tid1;
			if(i!=0) return i;
			i = pos1 - o.pos1;
			if(i!=0) return i;
			i = tid2 - o.tid2;
			if(i!=0) return i;
			i = pos2 - o.pos2;
			if(i!=0) return i;
			return 0;
			}
		
		int compare1(final Event o) {
			int i= compare2(o);
			if(i!=0) return i;
			return Long.compare(id, o.id);
			}
		}
	private static class EventCodec extends AbstractDataCodec<Event>
		{
		@Override
		public Event decode(final DataInputStream dis) throws IOException {
			final Event evt=new Event();
			try {
				evt.partition = dis.readUTF();
			} catch(IOException err) { return null;}
			evt.tid1=dis.readInt();
			evt.pos1=dis.readInt();
			evt.tid2=dis.readInt();
			evt.pos2=dis.readInt();		
			evt.id = dis.readLong();
			return evt;
			}
		@Override
		public void encode(final DataOutputStream dos, final Event o) throws IOException {
			dos.writeUTF(o.partition);
			dos.writeInt(o.tid1);
			dos.writeInt(o.pos1);
			dos.writeInt(o.tid2);
			dos.writeInt(o.pos2);
			dos.writeLong(o.id);
			}
		@Override
		public AbstractDataCodec<Event> clone() {
			return new EventCodec();
			}
		}
	
	private abstract class Report
		{
		protected final SAMSequenceDictionary dict;
		Report(final SAMSequenceDictionary dict) {
			this.dict=dict;
		}
		abstract void write(final List<Event> events) throws IOException,XMLStreamException;
		abstract void close() throws IOException,XMLStreamException;
		}
	private class TextReport extends Report
		{
		final PrintWriter w;
		TextReport(final File filename,final SAMSequenceDictionary dict) throws IOException {
			super(dict);
			this.w= (filename==null?new PrintWriter(stdout()):IOUtils.openFileForPrintWriter(filename));
			this.w.println("#"+ samRecordPartition.name() +"\tcontig1\tpos1\tcontig2\tpos2\tcount");
			}
		@Override
		void write(final List<Event> events) throws IOException,XMLStreamException
			{
			if(events.isEmpty()) return;
			final Event first=events.get(0);
			w.println(
				first.partition+"\t"+	
				dict.getSequence(first.tid1).getSequenceName()+"\t"+
				first.pos1+"\t"+
				dict.getSequence(first.tid2).getSequenceName()+"\t"+
				first.pos2+"\t"+
				events.size()
				);
			}
		@Override
		void close() throws IOException,XMLStreamException
			{
			w.flush();
			w.close();
			}
		}
	private class XMLReport extends Report
		{
		private String prevPartition = null;
		private final XMLStreamWriter w;
		XMLReport(File filename,final SAMSequenceDictionary dict) throws IOException,XMLStreamException {
			super(dict);
			final XMLOutputFactory xof =  XMLOutputFactory.newFactory();
			if(filename==null) {
				w =  xof.createXMLStreamWriter(stdout());
				}
			else
				{
				w =  xof.createXMLStreamWriter(new PrintWriter(filename,"UTF-8"));
				}
		
			
			w.writeStartDocument("UTF-8", "1.0");
			w.writeStartElement("translocations");
		
			long n=0;
			long genome_size = dict.getReferenceLength();
			w.writeStartElement("dictionary");
			this.w.writeAttribute("size",String.valueOf(dict.size()));
			this.w.writeAttribute("length",String.valueOf(genome_size));
			for(int i=0;i<dict.size();i++)
				{
				final SAMSequenceRecord rec=this.dict.getSequence(i);
				this.w.writeStartElement("contig");
				this.w.writeAttribute("name",rec.getSequenceName());
				this.w.writeAttribute("tid",String.valueOf(i));
				this.w.writeAttribute("length",String.valueOf(rec.getSequenceLength()));
				this.w.writeAttribute("index",String.valueOf(n));
				this.w.writeAttribute("index-fraction",String.valueOf(n/(double)genome_size));
				this.w.writeCharacters(rec.getSequenceName());
				this.w.writeEndElement();
				n+=rec.getSequenceLength();
				}
			w.writeEndElement();
			w.writeStartElement("partitions");
			}
		
		private void writeSplit(final String tag,int tid,int pos)  throws IOException,XMLStreamException {
			final SAMSequenceRecord rec=this.dict.getSequence(tid);
			this.w.writeStartElement(tag);
			this.w.writeAttribute("contig",rec.getSequenceName());
			this.w.writeAttribute("tid", String.valueOf(tid));
			this.w.writeAttribute("pos", String.valueOf(pos));
			long n=0;
			for(int i=0;i<tid;i++) n+=this.dict.getSequence(i).getSequenceLength();
			n+=pos;
			this.w.writeAttribute("index", String.valueOf(n));
			
			this.w.writeCharacters(rec.getSequenceName());
			this.w.writeEndElement();
			}
		
		void write(final List<Event> events) throws IOException,XMLStreamException
			{
			if(events.isEmpty()) return;
			final Event first=events.get(0);
			if(!first.partition.equals(prevPartition)) 
				{
				if(prevPartition!=null)
					{
					this.w.writeEndElement();
					}
				this.w.writeStartElement("partition");
				this.w.writeAttribute("name", first.partition);
				this.w.writeAttribute("partition-type", samRecordPartition.name());
				
				
				
				this.prevPartition = first.partition;
				}
			
			this.w.writeStartElement("event");
			this.w.writeAttribute("count", String.valueOf(events.size()));
			writeSplit("start",first.tid1,first.pos1);
			writeSplit("end",first.tid2,first.pos2);
			this.w.writeEndElement();
				
			}
		@Override
		void close() throws IOException, XMLStreamException {
			if(prevPartition!=null)
				{
				this.w.writeEndElement();
				}
			w.writeEndElement();//partitions
			w.writeEndElement();
			w.writeEndDocument();
			w.flush();
			w.close();
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.round_loc<=0) {
			LOG.error("round_loc <=0 ("+round_loc+")");
			return -1;
		}
		long id_generator=0L;
		ConcatSam.ConcatSamIterator samIter = null;
		SortingCollection<Event> events = null;
		CloseableIterator<Event> evtIter=null;
		try {
			IntFunction<Integer> roundPosition = POS->POS-POS%round_loc;
			
			samIter = new ConcatSam.Factory().setInterval(this.region_str).open(args);
			final SAMSequenceDictionary dict = samIter.getFileHeader().getSequenceDictionary();
			if(dict.size()<2) {
				LOG.error("Not enough contigs in sequence dictionary. Expected at least 2.");
				return -1;
			}
 			
			
			events = SortingCollection.newInstance(Event.class,
					new EventCodec(),
					(E1,E2)->E1.compare1(E2),
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			events.setDestructiveIteration(true);
			
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(dict).logger(LOG);
			while(samIter.hasNext()) {
				final SAMRecord rec = progress.watch(samIter.next());
				if(rec.getReadUnmappedFlag()) continue;
				if(!rec.getReadPairedFlag()) continue;
				if(rec.getMateUnmappedFlag()) continue;
				if(this.samRecordFilter.filterOut(rec)) continue;
				
				final int tid1 =  rec.getReferenceIndex();
				final int tid2 =  rec.getMateReferenceIndex();
				if(tid1==tid2) continue;
				
				final Event event = new Event();
				event.partition = this.samRecordPartition.getPartion(rec, "N/A");
				event.id = ++id_generator;
				if(tid1<tid2) {
					event.tid1 = tid1;
					event.pos1 = roundPosition.apply(rec.getStart());
					event.tid2 = tid2;
					event.pos2 = roundPosition.apply(rec.getMateAlignmentStart());
					}
				else
					{
					event.tid2 = tid1;
					event.pos2 = roundPosition.apply(rec.getStart());
					event.tid1 = tid2;
					event.pos1 = roundPosition.apply(rec.getMateAlignmentStart());
					}
				events.add(event);
				if(event.id%1000L==0) {
					LOG.info("number of translocation events : "+event.id);
					}
				}
			progress.finish();
			samIter.close();samIter=null;
			events.doneAdding();
			final Report report;
			if(this.xml_output)
				{
				report= new XMLReport(this.outputFile, dict);
				}
			else
				{
				report = new TextReport(this.outputFile,dict);
				}
			
			evtIter = events.iterator();
			final EqualRangeIterator<Event> eq=new EqualRangeIterator<>(evtIter, (E1,E2)->E1.compare2(E2));
			while(eq.hasNext())
				{
				final List<Event> eventList = eq.next();
				if(eventList.size()< this.min_number_of_event) continue;
				report.write(eventList);
				}
			report.close();
			eq.close();
			evtIter.close();evtIter=null;
			events.cleanup();
			events=null;
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(evtIter);
			if(events!=null) events.cleanup();
			events=null;
			CloserUtil.close(samIter);
			}
		}
	
	public static void main(String[] args) {
		new SamTranslocations().instanceMainWithExit(args);

	}

}
