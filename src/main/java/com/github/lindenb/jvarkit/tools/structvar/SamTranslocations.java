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
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
/**

BEGIN_DOC

## converting to SVG

with the following XSLT stylesheet (last updated : 2017-11-30 )

```xslt
<?xml version='1.0' encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl='http://www.w3.org/1999/XSL/Transform' 
	xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
	xmlns:svg="http://www.w3.org/2000/svg"
	xmlns="http://www.w3.org/2000/svg"
        xmlns:xlink="http://www.w3.org/1999/xlink"
        xmlns:str="http://exslt.org/strings"
        xmlns:math="http://exslt.org/math"
	extension-element-prefixes="str math"
	version='1.0'>
<xsl:output method="xml"  encoding="UTF-8"/>
<xsl:variable name="nsamples" select="count(translocations/partitions/partition)"/>
<xsl:variable name="cols" select="floor(math:sqrt($nsamples))"/>
<xsl:variable name="W" select="300"/>
<xsl:variable name="max-count">
  <xsl:for-each select="translocations/partitions/partition/event">
    <xsl:sort select="@count" data-type="number" order="descending" />
    <xsl:if test="position() = 1">
      <xsl:value-of select="number(@count)" />
    </xsl:if>
  </xsl:for-each>
</xsl:variable>



<xsl:template match="translocations">

<svg xmlns="http://www.w3.org/2000/svg">
  <xsl:attribute name="width">
  	<xsl:value-of select="$cols * $W"/>
  </xsl:attribute>
  <xsl:attribute name="height">
  	<xsl:value-of select="($cols +1) * $W"/>
  </xsl:attribute>
  <style>
        svg { fill:white;stroke:black;}
  	circle {fill:red; opacity:0.3; stroke:none;}
  	text {fill:gray;stroke:gray;font-size:10px;}
  	.c0 {fill:gray;opacity:0.05;}
  	.c1 {fill:white;opacity:0.05;}
  </style>
  <defs>
  	<g id="dict">
  		<xsl:apply-templates select="dictionary"/>
  	</g>
  </defs>
  <g>
  	<xsl:apply-templates select="partitions/partition"/>
  </g>
</svg>
</xsl:template>


<xsl:template match="dictionary">
<xsl:apply-templates select="contig"/>
<rect style="fill:none;stroke:gray;" x="0" y="0" width="{$W}"  height="{$W}">
</rect>
</xsl:template>

<xsl:template match="contig">
<xsl:variable name="clazz" select="concat('c',count(preceding-sibling::contig) mod 2)"/>
<xsl:variable name="genomelen" select="number(../@length)"/>
<xsl:variable name="x" select="(number(@index) div $genomelen) * $W"/>
<xsl:variable name="h" select="(number(@length) div $genomelen) * $W"/>

<rect class="{$clazz}"  x="0" y="{$x}" width="{$W}" height="{$h}">
</rect>
<rect class="{$clazz}"  y="0" x="{$x}" height="{$W}" width="{$h}">
</rect>
</xsl:template>


<xsl:template match="partition">
<xsl:variable name="idx" select="count(preceding-sibling::partition)"/>
<xsl:variable name="dx" select="($idx mod $cols) * $W"/>
<xsl:variable name="dy" select="floor(($idx div $cols)) * $W"/>
<g>
 <xsl:attribute name="transform">translate(<xsl:value-of select="$dx"/>,<xsl:value-of select="$dy"/>)</xsl:attribute>
 <title><xsl:value-of select="@name"/></title>
 <use x="0" y="0" href="#dict" />
 <text x="1" y="10"><xsl:value-of select="@name"/></text>
 <xsl:apply-templates select="event"/>
 
</g>
</xsl:template>

<xsl:template match="event">
<xsl:variable name="genomelen" select="number(../../../dictionary/@length)"/>
<xsl:variable name="cx" select="(number(start/@index) div $genomelen) * $W"/>
<xsl:variable name="cy" select="(number(end/@index) div $genomelen) * $W"/>
<xsl:variable name="r" select="0.1 + (number(@count) div $max-count) * 20.0"/>
<circle cx="{$cx}" cy="{$cy}" r="{$r}"/>
<circle cx="{$cy}" cy="{$cx}" r="{$r}"/>
</xsl:template>

</xsl:stylesheet>
```

END_DOC
*/
@Program(name="samtranslocations",
	description="Explore translocations between two chromosomes using discordant paired-end reads.",
	keywords={"sam","bam","xslt","xml"}
	)
public class SamTranslocations extends Launcher {
	private static final Logger LOG = Logger.build(SamTranslocations.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"--region","--interval"},description="Limit analysis to this interval. "+IntervalParser.OPT_DESC)
	private String region_str=null;
	@Parameter(names={"--filter"},description=SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamRecordJEXLFilter.StringConverter.class)
	private SamRecordFilter samRecordFilter = SamRecordJEXLFilter.buildDefault();
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
	    boolean clipped;
	    int count_partitions=0;
	    long id;
	    
		Event() {
			
			}
		
		int compareLoc(final Event o) {
			int i = tid1 - o.tid1;
			if(i!=0) return i;
			i = pos1 - o.pos1;
			if(i!=0) return i;
			i = tid2 - o.tid2;
			if(i!=0) return i;
			i = pos2 - o.pos2;
			if(i!=0) return i;
			return 0;
			}
		
		int comparePartitonLoc(final Event o) {
			final int i= partition.compareTo(o.partition);
			if(i!=0) return i;
			return compareLoc(o);
			}
		
		int comparePartitonLocId(final Event o) {
			final int i= comparePartitonLoc(o);
			if(i!=0) return i;
			return Long.compare(id, o.id);
			}
		int compareLocId(final Event o) {
			final int i= compareLoc(o);
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
			evt.clipped = dis.readBoolean();
			evt.count_partitions = dis.readInt();
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
			dos.writeBoolean(o.clipped);
			dos.writeInt(o.count_partitions);
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
			this.w.println("#"+ samRecordPartition.name() +"\tcontig1\tpos1\tcontig2\tpos2\tcount-events\tnumber_of_partition");
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
				events.size()+"\t"+
				first.count_partitions
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
		private int max_list_size = 0;
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
			this.w.writeAttribute("type", samRecordPartition.name());
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
			this.max_list_size = Math.max(max_list_size, events.size());
			final Event first=events.get(0);
			if(!first.partition.equals(prevPartition)) 
				{
				if(prevPartition!=null)
					{
					this.w.writeEndElement();
					}
				this.w.writeStartElement("partition");
				this.w.writeAttribute("name", first.partition);
				
				
				
				this.prevPartition = first.partition;
				}
			
			this.w.writeStartElement("event");
			this.w.writeAttribute("num-partitions",String.valueOf(first.count_partitions));
			this.w.writeAttribute("count", String.valueOf(events.size()));
			this.w.writeAttribute("count-clipped",String.valueOf( events.stream().filter(E->E.clipped).count()));
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
			
			this.w.writeStartElement("summary");
			this.w.writeEmptyElement("entry");
				this.w.writeAttribute("key", "max-count");
				this.w.writeAttribute("value",String.valueOf(this.max_list_size));
			this.w.writeEndElement();//summary
			
			w.writeEndElement();//translocations
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
		SortingCollection<Event> eventsLoc = null;
		SortingCollection<Event> eventsPartiton = null;
		CloseableIterator<Event> evtIter=null;
		try {
			IntFunction<Integer> roundPosition = POS->POS-POS%round_loc;
			
			samIter = new ConcatSam.Factory().setInterval(this.region_str).open(args);
			final SAMSequenceDictionary dict = samIter.getFileHeader().getSequenceDictionary();
			if(dict.size()<2) {
				LOG.error("Not enough contigs in sequence dictionary. Expected at least 2.");
				return -1;
			}
			
			eventsLoc = SortingCollection.newInstance(Event.class,
					new EventCodec(),
					(E1,E2)->E1.compareLocId(E2),
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			eventsLoc.setDestructiveIteration(true);
			
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
				event.clipped = rec.getCigar()!=null && rec.getCigar().isClipped();
				eventsLoc.add(event);
				if(event.id%1000L==0) {
					LOG.info("number of translocation events : "+event.id);
					}
				}
			progress.finish();
			samIter.close();samIter=null;
			eventsLoc.doneAdding();
			
			
			eventsPartiton = SortingCollection.newInstance(Event.class,
					new EventCodec(),
					(E1,E2)->E1.comparePartitonLocId(E2),
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			eventsPartiton.setDestructiveIteration(true);

			
			evtIter = eventsLoc.iterator();
			EqualRangeIterator<Event> eq=new EqualRangeIterator<>(evtIter, (E1,E2)->E1.compareLoc(E2));
			while(eq.hasNext())
				{
				final List<Event> eventList = eq.next();
				for(final Event evt:eventList) {
					evt.count_partitions = eventList.size();
					eventsPartiton.add(evt);
					}				
				}
			eq.close();
			eq=null;
			evtIter.close();evtIter=null;
			try { eventsLoc.cleanup();} catch(Throwable err){}
			eventsLoc=null;
			
			final Report report;
			if(this.xml_output)
				{
				report= new XMLReport(this.outputFile, dict);
				}
			else
				{
				report = new TextReport(this.outputFile,dict);
				}
			
			eventsPartiton.doneAdding();
			evtIter = eventsPartiton.iterator();
			eq=new EqualRangeIterator<>(evtIter, (E1,E2)->E1.comparePartitonLoc(E2));
			while(eq.hasNext())
				{
				final List<Event> eventList = eq.next();
				if(eventList.size()< this.min_number_of_event) continue;
				report.write(eventList);
				}
			report.close();
			eq.close();
			evtIter.close();evtIter=null;
			eventsPartiton.cleanup();
			eventsPartiton=null;
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(evtIter);
			if(eventsPartiton!=null) try {eventsPartiton.cleanup();} catch(Exception err){}
			eventsPartiton=null;
			CloserUtil.close(samIter);
			}
		}
	
	public static void main(String[] args) {
		new SamTranslocations().instanceMainWithExit(args);

	}

}
