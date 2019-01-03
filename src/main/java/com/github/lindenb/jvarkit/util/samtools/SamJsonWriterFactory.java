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
package com.github.lindenb.jvarkit.util.samtools;

import java.io.IOException;
import java.io.Writer;
import java.util.List;

import com.google.gson.stream.JsonWriter;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.ProgressLoggerInterface;
import htsjdk.samtools.util.RuntimeIOException;

public class SamJsonWriterFactory {
	private boolean printHeader = true;
	private boolean printReadName = true;
	private boolean printReadSequence = true;
	private boolean printReadQualities = true;
	private boolean printMate = true;
	private boolean printAttributes = true;
	private boolean expandflag = false;
	private boolean expandcigar = false;
	private boolean closeStream = true;
	
	
	public SamJsonWriterFactory printHeader(boolean v) { this.printHeader=v; return this;}
	public SamJsonWriterFactory printReadName(boolean v) { this.printReadName=v; return this;}
	public SamJsonWriterFactory printReadSequence(boolean v) { this.printReadSequence=v; return this;}
	public SamJsonWriterFactory printReadQualities(boolean v) { this.printReadQualities=v; return this;}
	public SamJsonWriterFactory printMate(boolean v) { this.printMate=v; return this;}
	public SamJsonWriterFactory printAttributes(boolean v) { this.printAttributes=v; return this;}
	public SamJsonWriterFactory expandFlag(boolean v) { this.expandflag=v; return this;}
	public SamJsonWriterFactory expandCigar(boolean v) { this.expandcigar=v; return this;}
	public SamJsonWriterFactory closeStreamAtEnd(boolean v) { this.closeStream=v; return this;}	
	
	public SamJsonWriterFactory() {
	}
	
	public static SamJsonWriterFactory newInstance() { return new SamJsonWriterFactory();}
	
	public SAMFileWriter open(final SAMFileHeader header,final Writer w) {
		return  open(header,new JsonWriter(w));
		}

	
	public SAMFileWriter open(final SAMFileHeader header,final JsonWriter w) {
		return new JSONWriter(header, w);
		}

private class JSONWriter implements SAMFileWriter
	{
	private final boolean printHeader = SamJsonWriterFactory.this.printHeader;
	private final boolean printReadName = SamJsonWriterFactory.this.printReadName;
	private final boolean printReadSequence = SamJsonWriterFactory.this.printReadSequence;
	private final boolean printReadQualities = SamJsonWriterFactory.this.printReadQualities;
	private final boolean printMate = SamJsonWriterFactory.this.printMate;
	private final boolean printAttributes = SamJsonWriterFactory.this.printAttributes;
	private final boolean expandflag = SamJsonWriterFactory.this.expandflag;
	private final boolean expandcigar = SamJsonWriterFactory.this.expandcigar;
	private final boolean closeStream = SamJsonWriterFactory.this.closeStream;
	private SAMFileHeader header = null;
	private ProgressLoggerInterface progress;
	private JsonWriter w;
	
	
	JSONWriter(final SAMFileHeader header,JsonWriter w) {
		this.w=w;
		this.header=header;
		
		
		
		try {
			if(this.printHeader ) {
				w.beginObject();
				w.name("header");
				w.beginObject();
				w.name("version");
				w.value(header.getVersion());
				w.name("sortorder");
				w.value(header.getSortOrder().name());
				w.name("dict");
				final SAMSequenceDictionary dict=header.getSequenceDictionary();
				if(dict==null)
					{
					w.nullValue();
					}
				else
					{
					w.beginArray();
					for(int i=0;i< dict.size();++i)
						{
						final SAMSequenceRecord rec=dict.getSequence(i);
						w.beginObject();
						w.name("name");
						w.value(rec.getSequenceName());
						w.name("length");
						w.value(rec.getSequenceLength());
						if(rec.getAssembly()!=null)
							{
							w.name("assembly");
							w.value(rec.getAssembly());
							}
						w.endObject();
						}
					w.endArray();
					}
				w.endObject();
				w.name("reads");
				}
			w.beginArray();
			crlf();
			
		} catch(IOException err) {
			throw new RuntimeIOException();
		}
	}
	
	private void crlf() throws IOException
		{
		//if(this.printcrlf) this.w.jsonValue("\n");
		}
	
	private void printObject(Object o) throws IOException
		{
		if(o==null)
			{
			w.nullValue();
			}
		else if(o instanceof List)
			{
			@SuppressWarnings("rawtypes")
			List<?> L=(List)o;
			w.beginArray();
			for(int i=0;i<L.size();++i)
				{
				printObject(L.get(i));
				}
			w.endArray();
			}
		else if(o.getClass().isArray())
			{
			w.beginArray();
			if(o instanceof byte[])
				{
				byte L[]=(byte[])o;
				for(int i=0;i< L.length;++i)
					{
					printObject((short)L[i]);
					}
				}
			else if(o instanceof short[])
				{
				short L[]=(short[])o;
				for(int i=0;i< L.length;++i)
					{
					printObject(L[i]);
					}
				}
			else if(o instanceof int[])
				{
				 int L[]=( int[])o;
				for(int i=0;i< L.length;++i)
					{
					printObject(L[i]);
					}
				}
			else if(o instanceof long[])
				{
				long L[]=( long[])o;
				for(int i=0;i< L.length;++i)
					{
					printObject(L[i]);
					}
				}
			else if(o instanceof float[])
				{
				float L[]=( float[])o;
				for(int i=0;i< L.length;++i)
					{
					printObject(L[i]);
					}
				}
			else if(o instanceof double[])
				{
				double L[]=( double[])o;
				for(int i=0;i< L.length;++i)
					{
					printObject(L[i]);
					}
				}
			else
				{
				Object L[]=( Object[])o;
				for(int i=0;i< L.length;++i)
					{
					printObject(L[i]);
					}
				}
			w.endArray();
			}
		else if(o instanceof Number )
			{
			w.value(Number.class.cast(o));
			}
		else if(o instanceof Boolean)
			{
			w.value(Boolean.class.cast(o));
			}
		else
			{
			w.value(String.valueOf(o));
			}
		}

	
	@Override
	public void addAlignment(final SAMRecord rec) {
		if(this.progress!=null) this.progress.record(rec);
		try {
			this.w.beginObject();
			if(this.printReadName) {
				this.w.name("name");
				this.w.value(rec.getReadName());
			}
			this.w.name("flag");
			if(this.expandflag)
				{
				w.beginObject();
				for(SAMFlag flg: SAMFlag.values()) 
					{
					w.name(flg.name());
					w.value(flg.isSet(rec.getFlags()));
					}
				w.endObject();
				}
			else
				{
				this.w.value(rec.getFlags());
				}
			
			if(rec.getReferenceName()!=null)
				{
				this.w.name("ref");
				this.w.value(rec.getContig());
				}
			
			this.w.name("pos");
			this.w.value(rec.getAlignmentStart());
			
			if(!rec.getReadUnmappedFlag()) {
				this.w.name("mapq");
				this.w.value(rec.getMappingQuality());
				
				this.w.name("cigar");

				if(this.expandcigar) {
					final Cigar  cigar= rec.getCigar();
					if(cigar==null) {
						this.w.nullValue();
					} else
					{
						this.w.beginArray();
						for(final CigarElement ce: cigar) {
							this.w.beginObject();
							this.w.name("op");
							this.w.value(ce.getOperator().name());
							this.w.name("len");
							this.w.value(ce.getLength());
							this.w.endObject();
						}
						
						this.w.endArray();
					}
				} else
				{
					this.w.value(rec.getCigarString());
				}
				
			}

			
			
			
			if(this.printMate && rec.getReadPairedFlag()){
				this.w.name("len");
				this.w.value(rec.getInferredInsertSize());
				
				if(rec.getMateReferenceName()!=null)
					{
					this.w.name("materef");
					this.w.value(rec.getMateReferenceName());
					this.w.name("matepos");
					this.w.value(rec.getMateAlignmentStart());
					}
				}
			
				
			if(this.printReadSequence) {
				this.w.name("sequence");
				this.w.value(rec.getReadString());
			}
			if(this.printReadQualities) {
				this.w.name("qualities");
				this.w.value(rec.getBaseQualityString());
			}
			
			if(this.printAttributes) {
				this.w.name("atts");
				this.w.beginArray();
				List<SAMRecord.SAMTagAndValue> atts=rec.getAttributes();
				for(int i=0;i<atts.size();++i)
					{
					
					final SAMRecord.SAMTagAndValue att=atts.get(i);
					this.w.beginObject();
					this.w.name("name");
					this.w.value(att.tag);
					this.w.name("value");
					printObject(att.value);
					this.w.endObject();
					}
				this.w.endArray();
				}
			
			
			this.w.endObject();
			crlf();
		} catch(IOException err) {
			throw new RuntimeIOException();
		}
		
	}

	@Override
	public SAMFileHeader getFileHeader() {
		return header;
	}

	@Override
	public void setProgressLogger(final ProgressLoggerInterface progress) {
		this.progress=progress;
		
	}

	@Override
	public void close() {
		try {
			w.endArray();
			if(this.printHeader ) {
				w.endObject();
			}
			crlf();
			w.flush();
		} catch(IOException err) {
			throw new RuntimeIOException();
		}
		if(this.closeStream) CloserUtil.close(this.w);
		}
	
	}
}
