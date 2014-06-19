package com.github.lindenb.jvarkit.util.picard;

import java.io.PrintWriter;
import java.util.List;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.ProgressLoggerInterface;

public class SamJsonWriter implements SAMFileWriter
	{
	private PrintWriter out;
	private SAMFileHeader header;
	private long count=0L;
	private boolean crlf=false;
	private boolean header_printed=false;
	private boolean print_header=true;
	@SuppressWarnings("unused")
	private ProgressLoggerInterface progress;
	public SamJsonWriter(PrintWriter out,SAMFileHeader header)
		{
		this.out=out;
		this.header=header;
		}
	
	private void printHeader()
		{
		if(header_printed) return;
		header_printed=true;
		if(print_header)
			{
			out.print("{\"version\":\"");
			out.print(getFileHeader().getVersion());
			out.print("\",\"sortorder\":\"");
			out.print(getFileHeader().getSortOrder());
			out.print("\",\"dict\":");
			
			SAMSequenceDictionary dict=getFileHeader().getSequenceDictionary();
			if(dict==null)
				{
				out.print("null");
				}
			else
				{
				out.print("[");
				for(int i=0;i< dict.size();++i)
					{
					SAMSequenceRecord rec=dict.getSequence(i);
					if(i>0) out.print(',');
					out.print("{\"name\":\"");
					out.print(rec.getSequenceName());
					out.print("\",\"length\":");
					out.print(rec.getSequenceLength());
					out.print("}");
					}
				
				out.print("]");
				}
			out.print(",\"alignments\":");
			}
		out.print("[");
		}
	public void setAddCarriageReturn(boolean crlf)
		{		
		this.crlf=crlf;
		}
	
	public void setPrintHeader(boolean print_header)
		{
		this.print_header = print_header;
		}
	
	@Override
	public void addAlignment(SAMRecord rec)
		{
		if(out==null) return;//out closed
		printHeader();
		if(count>0L) out.print(",");
		if(crlf) out.println();
		out.print("{\"name\":\"");
		out.print(rec.getReadName());
		out.print("\",\"flag\":");
		out.print(rec.getFlags());
		out.print(",\"ref\":\"");
		out.print(rec.getReferenceName());
		out.print("\",\"start\":");
		out.print(rec.getAlignmentStart());
		out.print(",\"mapq\":");
		out.print(rec.getMappingQuality());
		out.print(",\"cigar\":");
		if(rec.getReadUnmappedFlag())
			{
			out.print("null");
			}
		else
			{
			out.print("\"");
			out.print(rec.getCigarString());
			out.print("\"");
			}
		if(rec.getReadPairedFlag())
			{
			out.print(",\"materef\":\"");
			out.print(rec.getMateReferenceName());
			out.print("\",\"matestart\":");
			out.print(rec.getMateAlignmentStart());
			
			if(rec.getInferredInsertSize()!=0)
				{
				out.print(",\"len\":");
				out.print(rec.getInferredInsertSize());
				}
			}
		out.print(",\"seq\":\"");
		out.print(rec.getReadString());
		out.print("\",\"qual\":\"");
		out.print(rec.getBaseQualityString());
		out.print("\",\"atts\":[");
		List<SAMRecord.SAMTagAndValue> atts=rec.getAttributes();
		for(int i=0;i<atts.size();++i)
			{
			if(i>0) out.print(",");
			SAMRecord.SAMTagAndValue att=atts.get(i);
			out.print("{\"name\":\"");
			out.print(att.tag);
			out.print("\",\"value\":");
			printObject(att.value);
			out.print("}");
			}
		out.print("]");
		
		out.print("}");
		++count;
		}

	private void printObject(Object o)
		{
		if(o==null)
			{
			out.print("null");
			}
		else if(o instanceof List)
			{
			@SuppressWarnings("rawtypes")
			List<?> L=(List)o;
			out.print("[");
			for(int i=0;i<L.size();++i)
				{
				if(i>0) out.print(",");
				printObject(L.get(i));
				}
			out.print("]");
			}
		else if(o.getClass().isArray())
			{
			out.print("[");
			if(o instanceof byte[])
				{
				byte L[]=(byte[])o;
				for(int i=0;i< L.length;++i)
					{
					if(i>0) out.print(",");
					printObject((short)L[i]);
					}
				}
			else if(o instanceof short[])
				{
				short L[]=(short[])o;
				for(int i=0;i< L.length;++i)
					{
					if(i>0) out.print(",");
					printObject(L[i]);
					}
				}
			else if(o instanceof int[])
				{
				 int L[]=( int[])o;
				for(int i=0;i< L.length;++i)
					{
					if(i>0) out.print(",");
					printObject(L[i]);
					}
				}
			else if(o instanceof long[])
				{
				long L[]=( long[])o;
				for(int i=0;i< L.length;++i)
					{
					if(i>0) out.print(",");
					printObject(L[i]);
					}
				}
			else if(o instanceof float[])
				{
				float L[]=( float[])o;
				for(int i=0;i< L.length;++i)
					{
					if(i>0) out.print(",");
					printObject(L[i]);
					}
				}
			else if(o instanceof double[])
				{
				double L[]=( double[])o;
				for(int i=0;i< L.length;++i)
					{
					if(i>0) out.print(",");
					printObject(L[i]);
					}
				}
			else
				{
				Object L[]=( Object[])o;
				for(int i=0;i< L.length;++i)
					{
					if(i>0) out.print(",");
					printObject(L[i]);
					}
				}
			out.print("]");
			}
		else if(o instanceof Number ||o instanceof Boolean)
			{
			out.print(String.valueOf(o));
			}
		else
			{
			out.print("\"");
			out.print(String.valueOf(o));
			out.print("\"");
			}
		}
	
	@Override
	public SAMFileHeader getFileHeader() {
		return this.header;
		}

	
	@Override
	public void close()
		{
		if(out==null) return;
		printHeader();
		if(crlf) out.println();
		out.print("]");
		if(print_header) out.print("}");
		out.flush();
		out=null;
		}

	@Override
	public void setProgressLogger(ProgressLoggerInterface progress) {
		this.progress=progress;
	}

	}
