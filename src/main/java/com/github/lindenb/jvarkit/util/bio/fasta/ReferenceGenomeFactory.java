/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util.bio.fasta;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLEncoder;
import java.util.OptionalDouble;
import java.util.OptionalInt;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.beust.jcommander.IStringConverter;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;

public class ReferenceGenomeFactory
implements IStringConverter<ReferenceGenome>  {
private static final Logger LOG = Logger.build(ReferenceGenomeFactory.class).make();

public static final String OPT_DESCRIPTION_FILE_ONLY="Indexed Genome Reference. "+
		"A fasta file that must be indexed with samtools faidx and with picard CreateSequenceDictionary."
		;

	
public static final String OPT_DESCRIPTION="Indexed Genome Reference. "+
			"It can be a the path to fasta file that must be indexed with samtools faidx and with picard CreateSequenceDictionary."
			+ " It can also be a BioDAS dsn url like `http://genome.cse.ucsc.edu/cgi-bin/das/hg19/` . BioDAS references are slower, but allow to work without a local reference file.";

/** jcommander stuff */
@Override
public ReferenceGenome convert(final String opt) {
		try {
			return this.open(opt);
			}
		catch(final IOException err) {
			LOG.error(err);
			throw new RuntimeIOException(err);
			}
	}

private final int DEFAULT_HALF_BUFFER_CAPACITY = 1_000_000;
private int half_buffer_capacity=DEFAULT_HALF_BUFFER_CAPACITY;
private boolean return_N_on_indexOutOfRange = false;
private boolean debug = false;
private boolean throwOnContigNotFound = false;
private boolean neverReturnNullContig = false;
private boolean disableDefaultAliase = false;

/** never return a null contig if it's not in the dict, instead return a 0-length contig that will always return 'N' for 'charAt(idx)' */
public ReferenceGenomeFactory setNeverReturnNullContig(boolean neverReturnNullContig) {
	this.neverReturnNullContig = neverReturnNullContig;
	return this;
	}

public boolean isNeverReturnNullContig() {
	return neverReturnNullContig;
	}

/** throw an exception if contig is not found */
public ReferenceGenomeFactory setThrowOnContigNotFound(boolean throwOnContigNotFound) {
	this.throwOnContigNotFound = throwOnContigNotFound;
	return this;
	}

public boolean isThrowOnContigNotFound() {
	return throwOnContigNotFound;
	}

public ReferenceGenomeFactory setBufferSize(int size) {
	this.half_buffer_capacity =Math.max(1,size);
	return this;
}

public int getBufferSize() {
	return this.half_buffer_capacity;
	}

public ReferenceGenomeFactory setReturnBaseNOnIndexOutOfRange(boolean b) {
	this.return_N_on_indexOutOfRange = b;
	return this;
}

public boolean isReturnBaseNOnIndexOutOfRange() {
	return return_N_on_indexOutOfRange;
}

public ReferenceGenomeFactory setDebug(boolean b) {
	this.debug = b;
	return this;
}

public boolean isDebug() {
	return debug;
}
public ReferenceGenomeFactory setDisableDefaultAliases(boolean b) {
	this.disableDefaultAliase = b;
	return this;
}

public boolean isDisableDefaultAliases() {
	return disableDefaultAliase;
}



private class NullReferenceContig
	extends AbstractCharSequence
	implements ReferenceContig
	{
	private final SAMSequenceRecord ssr;
	NullReferenceContig(final String name) {
		this.ssr = new SAMSequenceRecord(name, 0);
		}
	
	@Override
	public SAMSequenceRecord getSAMSequenceRecord() {
		return this.ssr;
		}
	@Override
	public char charAt(int index) {
		return 'N';
		}
	@Override
	public GCPercent getGCPercent(final int start, final int end) {
		return new GCPercent() {
			@Override
			public int getATCount() {
				return 0;
				}
			@Override
			public OptionalInt getGCPercentAsInteger() {
				return OptionalInt.empty();
				}
			@Override
			public boolean isEmpty() {
				return true;
				}
			@Override
			public int getGCCount() {
				return 0;
				}
			@Override
			public OptionalDouble getGCPercent() {
				return OptionalDouble.empty();
				}
			@Override
			public int getStart() {
				return start+1;
				}
			@Override
			public int getEnd() {
				return end;
				}
			@Override
			public String getContig() {
				return NullReferenceContig.this.getContig();
				}
			@Override
			public int getAllCount() {
				return 1+end-start;
				}
			};
		}
	}


/** base class for a ReferenceGenome */
private abstract class AbstractReferenceGenome
	implements ReferenceGenome
	{
	protected SAMSequenceDictionary dictionary =null ;
	private ReferenceContig last_contig = null;

	protected abstract ReferenceContig create(final SAMSequenceRecord ssr);
	
	@Override
	public final ReferenceContig getContig(final String contigName) {
		if(this.last_contig!=null && last_contig.hasName(contigName)) {
			return this.last_contig;
			}
		final SAMSequenceRecord ssr=getDictionary().getSequence(contigName);
		if(ssr==null) {
			if(isNeverReturnNullContig()) {
				if(isDebug()) LOG.warn("returning non-null contig for "+contigName);
				this.last_contig =  new NullReferenceContig(contigName);
				return this.last_contig;
			}
			if(isThrowOnContigNotFound()) {
				throw new JvarkitException.ContigNotFoundInDictionary(contigName, getDictionary());
				}
			return null;
		}
		this.last_contig= create(ssr);
		return this.last_contig;
		}

	@Override
	public SAMSequenceDictionary getDictionary() {
		return this.dictionary;
		}
	@Override
	public int hashCode() {
		return getSource().hashCode();
		}
	
	@Override
	public boolean equals(final Object obj) {
		return obj==this;
		}
	
	@Override
	public String toString() {
		return getSource();
		}
	}


private static class GCPercentImpl implements ReferenceContig.GCPercent
	{
	final String contig;
	final int start1;
	final int end1;
	int count=0;
	int count_gc=0;
	int count_at=0;
	
	GCPercentImpl(String contig,int s1,int e1) {
		this.contig = contig;
		this.start1=s1;
		this.end1=e1;
		}
	@Override public int getAllCount() { return this.count;}
	@Override public int getGCCount() { return this.count_gc;}
	@Override public int getATCount(){ return this.count_at;}
	@Override
	public boolean isEmpty() { return this.count == 0; }
	@Override
	public OptionalDouble getGCPercent() {
		return ( this.count==0? 
				OptionalDouble.empty():
				OptionalDouble.of(this.count_gc/(double)this.count))
				;
		}
	@Override
	public OptionalInt getGCPercentAsInteger() {
		return ( this.count==0?
				OptionalInt.empty():
				OptionalInt.of((int)(getGCPercent().getAsDouble()*100.0)));
		}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + count;
		result = prime * result + count_at;
		result = prime * result + count_gc;
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (obj == null || !(obj instanceof GCPercentImpl)) {
			return false;
			}
		final GCPercentImpl other = (GCPercentImpl) obj;
		return	this.count==other.count &&
				this.count_at==other.count_at &&
				this.count_gc==other.count_gc;
		}
	@Override
	public String toString() {
		return "gc_percent("+this.count_gc+"/"+this.count+")="+this.getGCPercent();
		}
	@Override
	public String getContig() {
		return this.contig;
		}
	@Override
	public int getStart() {
		return this.start1;
		}
	@Override
	public int getEnd() {
		return this.end1;
		}
}


private abstract class AbstractReferenceContigImpl
	extends AbstractCharSequence
	implements ReferenceContig
	{
	private final ReferenceGenome owner;
	private final SAMSequenceRecord samSequenceRecord;
	private byte buffer[]=null;
	private int buffer_pos=-1;
	private final int half_buffer_capacity=ReferenceGenomeFactory.this.getBufferSize();
	
	protected abstract  byte[] refill(int start0,int end0);
	
	protected AbstractReferenceContigImpl(final ReferenceGenome owner,final SAMSequenceRecord ssr) {
		this.owner=owner;
		this.samSequenceRecord = ssr;
	}
	
	@Override
	public boolean hasName(final String name) {
		if(this.getContig().equals(name)) return true;
		final SAMSequenceRecord ssr2 = this.owner.getDictionary().getSequence(name);
		if(ssr2==null) return false;
		if(ssr2==this.samSequenceRecord) return true;
		return ssr2.getSequenceName().equals(this.getContig());
		}

	@Override
	public final SAMSequenceRecord getSAMSequenceRecord() {
		return this.samSequenceRecord;
		}
	
	@Override
	public final char charAt(int index0) {
		if(index0<0 || index0 >= length())
			{
			if(ReferenceGenomeFactory.this.isReturnBaseNOnIndexOutOfRange()) {
				if(isDebug()) LOG.debug("index out of range "+index0);
				return 'N';
				}
			throw new IndexOutOfBoundsException("index:"+index0);
			}
		if(buffer!=null && index0>=buffer_pos && index0-buffer_pos < buffer.length)
			{
			return (char)buffer[index0-buffer_pos];
			}
		final int minStart=Math.max(0, index0-half_buffer_capacity);
		final int maxEnd=Math.min(minStart+2*half_buffer_capacity,this.length());
		if(isDebug()) {
			LOG.debug("Refill "+minStart+" to "+maxEnd);
			}
		this.buffer=refill(minStart,maxEnd);
		this.buffer_pos=minStart;
		return (char)buffer[index0-minStart];
		}
	
	@Override
	public GCPercent getGCPercent(final int start,final int end) {
		final int L=this.length();
		final GCPercentImpl gcp = new GCPercentImpl(
				this.getContig(),
				start+1,
				Math.min(end, L)
				);
		for(int i=start;i< end && i< L;++i) {
			gcp.count++;
			switch(this.charAt(i)) {
				case 'c': case 'C':
				case 'g': case 'G':
				case 's': case 'S':gcp.count_gc++; break;
				case 'a': case 'A':
				case 't': case 'T':
				case 'w': case 'W':gcp.count_at++; break;
				}
			}
		return gcp;
		}
	}
	
private  class ReferenceGenomeImpl
	extends AbstractReferenceGenome
	{
	private  class ReferenceContigImpl
		extends AbstractReferenceContigImpl
		{
		final ReferenceSequence referenceSequence;

		ReferenceContigImpl(final SAMSequenceRecord ssr) {
			super(ReferenceGenomeImpl.this,ssr);
			this.referenceSequence  = ReferenceGenomeImpl.this.indexedFastaSequenceFile.getSequence(ssr.getSequenceName());
			if(this.referenceSequence==null) throw new IllegalStateException();
			}

		@Override
		protected byte[] refill(int minStart,int maxEnd) {
			return ReferenceGenomeImpl.this.indexedFastaSequenceFile.getSubsequenceAt(
					getContig(),
					minStart+1,
					Math.min(maxEnd,this.length())
					).getBases();
			}
		
		}
	
	private final File fastaFile;
	private IndexedFastaSequenceFile indexedFastaSequenceFile;
	ReferenceGenomeImpl(final File fastaFile) throws IOException
		{
		this.fastaFile = fastaFile;
		IOUtil.assertFileIsReadable(fastaFile);
		this.indexedFastaSequenceFile = new IndexedFastaSequenceFile(fastaFile);
		super.dictionary = this.indexedFastaSequenceFile.getSequenceDictionary();
		if(super.dictionary==null) {
			throw new JvarkitException.FastaDictionaryMissing(fastaFile);
			}
		if(!ReferenceGenomeFactory.this.isDisableDefaultAliases()) {
			ContigNameConverter.setDefaultAliases(super.dictionary);
			}
		}
	@Override
	public String getSource() {
		return this.fastaFile.toString();
		}
	@Override
	protected ReferenceContig create(SAMSequenceRecord ssr) {
		return new ReferenceContigImpl(ssr);
		}
	
	@Override
	public void close() throws IOException {
		CloserUtil.close(this.indexedFastaSequenceFile);
		}
	
	}




private class DasGenomeImpl extends AbstractReferenceGenome
	{
	final String basedasurl;
	final XMLInputFactory xmlInputFactory =XMLInputFactory.newFactory();
	
	private class DasContig extends AbstractReferenceContigImpl
		{
		DasContig(final SAMSequenceRecord ssr) {
			super(DasGenomeImpl.this,ssr);
			}
				
		@Override
		protected byte[] refill(int start0, int end0) {
			InputStream in = null;
			XMLEventReader xef = null;
			try {
				final String dna_url = DasGenomeImpl.this.basedasurl+"/dna?segment="+
						URLEncoder.encode(getContig(), "UTF-8")+","+
						(start0+1)+","+end0;
				if(isDebug()) LOG.debug(dna_url);
				final QName DNA= new QName("DNA");
				in = new URL(dna_url).openStream();
				xef = DasGenomeImpl.this.xmlInputFactory.createXMLEventReader(in);
				while(xef.hasNext())
					{
					XMLEvent evt=xef.nextEvent();
					if(evt.isStartElement() && evt.asStartElement().getName().equals(DNA)) 
						{
						ByteArrayOutputStream baos = new ByteArrayOutputStream(Math.max(1,end0-start0));
						while(xef.hasNext())
							{
							evt=xef.nextEvent();
							if(evt.isCharacters()) {
								final String sequence = evt.asCharacters().getData();
								for(int i=0;i< sequence.length();i++)
									{
									if(Character.isWhitespace(sequence.charAt(i))) continue;
									baos.write((byte)sequence.charAt(i));
									}
								
								}
							else if(evt.isEndElement())
								{
								baos.close();
								return baos.toByteArray();
								}
							else
								{
								throw new XMLStreamException(dna_url+ " : illegal state",evt.getLocation());
								}
							}
						throw new XMLStreamException(dna_url+ " : illegal state",evt.getLocation());
						}
					}
				throw new IllegalStateException(dna_url+ " : No <DNA> found");
				}
			catch(final Exception err)
				{
				LOG.error(err);
				throw new RuntimeIOException(err);
				}
			finally
				{
				CloserUtil.close(xef);
				CloserUtil.close(in);
				}
			}
		}
	
	
	DasGenomeImpl(String base) throws IOException {
		if(!base.endsWith("/")) base+="/";
		this.basedasurl = base;
		super.dictionary=new SAMSequenceDictionary();
		InputStream in = null;
		XMLEventReader xef = null;
		try {
			final QName SEGMENT= new QName("SEGMENT");
			final QName ATT_ID= new QName("id");
			final QName ATT_STOP= new QName("stop");
			final QName ATT_END= new QName("end");
			final QName ATT_LENGTH= new QName("length");
			final String entry_points_url = this.basedasurl+"entry_points"; 
			if(isDebug()) LOG.debug("parsing "+entry_points_url);
			in = new URL(entry_points_url).openStream();
			xef = this.xmlInputFactory.createXMLEventReader(in);
			while(xef.hasNext())
				{
				final XMLEvent evt=xef.nextEvent();
				if(!evt.isStartElement()) continue;
				final StartElement SE=evt.asStartElement();
				if(!SE.getName().equals(SEGMENT)) continue;
				Attribute att= SE.getAttributeByName(ATT_ID);
				if(att==null) throw new XMLStreamException(entry_points_url+":cannot get @id", SE.getLocation());
				final String id = att.getValue();
				att= SE.getAttributeByName(ATT_STOP);
				if(att==null) att=SE.getAttributeByName(ATT_END);
				if(att==null) att=SE.getAttributeByName(ATT_LENGTH);
				if(att==null) throw new XMLStreamException(entry_points_url+":cannot get @stop / @length / @end", SE.getLocation());
				final int length = Integer.parseInt(att.getValue());
				if(length<=0) throw new XMLStreamException("bad end "+length, SE.getLocation());
				super.dictionary.addSequence(new SAMSequenceRecord(id, length));
				}
			if(!ReferenceGenomeFactory.this.isDisableDefaultAliases()) {
				ContigNameConverter.setDefaultAliases(this.dictionary);
				}
			if(isDebug()) LOG.debug("dict in "+entry_points_url+" size : "+this.dictionary.size());
			}
		catch(final XMLStreamException err)
			{
			LOG.error(err);
			throw new IOException(err);
			}
		finally
			{
			CloserUtil.close(xef);
			CloserUtil.close(in);
			}
		}
	
	@Override
	public String getSource() {
		return this.basedasurl;
		}
	
	@Override
	protected ReferenceContig create(final SAMSequenceRecord ssr) {
		return new DasContig(ssr);
		}
	
	@Override
	public void close() throws IOException {
		
		}
	}

/** open any kind of reference */
public ReferenceGenome open(final String ref) throws IOException
	{
	if(StringUtil.isBlank(ref)) throw new IllegalArgumentException("null/empty arg");
	return (IOUtil.isUrl(ref))?
		openDAS(new URL(ref)):
		openFastaFile(new File(ref))
		;
	}
/** open a FASTA reference */
public ReferenceGenome openFastaFile(final File fastaFile) throws IOException
	{
	return new ReferenceGenomeImpl(fastaFile);
	}

/** open a DAS URL */
public ReferenceGenome openDAS(final URL dasUrl) throws IOException
	{
	return new DasGenomeImpl(dasUrl.toString());
	}
}
