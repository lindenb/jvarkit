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
package com.github.lindenb.jvarkit.tools.bam4deseq;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StreamTokenizer;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.bin.SamSequenceRecordBinMap;
import com.github.lindenb.jvarkit.util.command.Command;

public class HtSeqCount extends AbstractHtSeqCount
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(HtSeqCount.class);

	@Override
	public Command createCommand()
		{
		return new MyCommand();
		}
	
	static private class MyCommand extends AbstractHtSeqCount.AbstractHtSeqCountCommand
		{
	/** features we want in the GTF */
	private Set<String> _features=null;
	/** first Dict found, to print chrom names and compare with others */
	private SAMSequenceDictionary firstDict=null;
	/** mapping position to transcript */
	private SamSequenceRecordBinMap<Transcript> pos2transcript=null;
	/** all filenames */
	private List<String> filenames=new ArrayList<String>();
	
	/* current file scanner in filenames */
	private int file_index=0;
	/* user feature  in the GTF */
	/* all transcripts to print at the end */
	private HashMap<String, Transcript> name2transcript=new HashMap<String, Transcript>();
	/* first ouput line is header */
	private boolean print_header=false;

	
	private static class Transcript
		{
		String name;
		int count[];
		boolean bad_flag=false;//multiple chroms
		int tid;
		int start=Integer.MAX_VALUE;
		int end=0;
		}
	
	
	
	private void parseGTF(BufferedReader in,SAMSequenceDictionary dict) throws IOException
		{
		this.pos2transcript=new SamSequenceRecordBinMap<Transcript>(dict);
		final Pattern tab=Pattern.compile("[\t]");
		String line;
		while((line=in.readLine())!=null)
			{
			if(line.isEmpty() || line.startsWith("#")) continue;
			String tokens[]=tab.split(line);
			if(tokens.length<=8)
				{
				throw new IOException("Bad GTF line " +line);
				}
			if(!this._features.contains(tokens[2])) continue;
			String chrom=tokens[0];
			int tid=dict.getSequenceIndex(chrom);
			if(tid<0)
				{
				LOG.warn("Unknown chromosome "+chrom+" Ignoring "+line);
				continue;
				}
			String transcript_name=null;
			String info=tokens[8];
			StreamTokenizer scanner=new StreamTokenizer(new StringReader(info));
			scanner.wordChars('_', '_');
			scanner.wordChars('0', '9');
			while(scanner.nextToken()!=StreamTokenizer.TT_EOF)
				{
				String key=scanner.sval;
				//info("K "+key);
				if(scanner.nextToken()==StreamTokenizer.TT_EOF)
					{
					throw new IOException("Bad gtf in "+line+" after "+key);
					}
			
				String value=scanner.sval;
				//info("v "+value);
				if(key.equals(this.transcript_id))
					{
					transcript_name=value;
					break;
					}
				if(scanner.nextToken()==StreamTokenizer.TT_EOF) break;
				if(scanner.ttype!=';')
					{
					throw new IOException("Bad gtf in "+line+" expect ; ");
					}
				}
			if(transcript_name==null)
				{
				LOG.info("No "+transcript_id+" in "+line);
				continue;
				}
			
			
			Transcript transcript=name2transcript.get(transcript_name);
			if(transcript==null)
				{
				transcript=new Transcript();
				transcript.tid=tid;
				transcript.name=transcript_name;
				transcript.count=new int[filenames.size()];
				Arrays.fill(transcript.count,0);
				name2transcript.put(transcript_name,transcript);
				if(name2transcript.size()%1000==0)
					{
					LOG.info("Transcripts: "+name2transcript.size()+" "+transcript_name);
					}
				}
			else
				{
				if(transcript.tid!=tid)
					{
					LOG.warn("Multiple chromosomes for "+transcript);
					transcript.bad_flag=true;
					continue;
					}
				}
			
			int start1=Integer.parseInt(tokens[3]);
			int end1=Integer.parseInt(tokens[4]);
			transcript.start=Math.min(transcript.start,start1);
			transcript.end=Math.max(transcript.end,end1);
			
			
			this.pos2transcript.put(tid, start1-1, end1,transcript);
			
			/**
			boolean ok=false;
			Iterator<Transcript>x=this.pos2transcript.overlapping(tid, start1-1, end1);
			while(x.hasNext()) if(transcript==x.next()) {ok=true;}
			if(!ok) throw new IllegalStateException("boum");*/
			}
		LOG.info("Done Reading transcripts:"+name2transcript.size());
		}
	
	private void touch(int tid,int start1,int end1)
		{
		Set<String> seen=null;
		Iterator<Transcript> iter=this.pos2transcript.overlapping(tid, start1-1, end1);
		while(iter.hasNext())
			{
			Transcript transcript=iter.next();
			if(seen==null) seen=new HashSet<String>();
			if(!seen.add(transcript.name)) continue;
			transcript.count[this.file_index]++;
			}
		}
	
	private void run(SamReader sfr) throws IOException
		{
		SAMFileHeader header=sfr.getFileHeader();
		
		if(this.file_index==0)
			{
			firstDict=header.getSequenceDictionary();
			
			LOG.info("Reading "+this.gtfFile);
			BufferedReader gtfIn=IOUtils.openFileForBufferedReading(this.gtfFile);
			parseGTF(gtfIn,header.getSequenceDictionary());
			gtfIn.close();
			}
		else
			{
			if(!SequenceUtil.areSequenceDictionariesEqual(firstDict, header.getSequenceDictionary()))
				{
				throw new IOException("not same sequence dictionaires between "+
					filenames.get(0)+" && "+
					filenames.get(this.file_index)
					);
				}
			}
		long nReads=0L;
		SAMRecordIterator iter=sfr.iterator();
		while(iter.hasNext())
			{
			SAMRecord rec=iter.next();
			if(rec.getReadUnmappedFlag()) continue;
			if(rec.getDuplicateReadFlag()) continue;
			if(rec.getNotPrimaryAlignmentFlag()) continue;
			if(rec.getReadFailsVendorQualityCheckFlag()) continue;
			if(nReads++%1E6==0)
				{
				LOG.info("Read "+nReads+" in "+this.filenames.get(this.file_index));
				}
			touch(rec.getReferenceIndex(),
					rec.getAlignmentStart(),
					rec.getAlignmentEnd()
					);
			}
		iter.close();
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		throw new IllegalStateException("should never be called");
		}
	
	@Override
	public Collection<Throwable> call() throws Exception
			{
			if(super.userFeatures==null)
				{
				this._features=new HashSet<String>();
				for(String F:new String[]{"CDS","exon"})
					{
					LOG.info("Adding "+F+" as default feature.");
					this._features.add(F);
					}
				}
			else
				{
				this._features=new HashSet<String>(super.userFeatures);
				}
			if( this.gtfFile==null)
				{
				return wrapException("undefined GTF file");
				}
			final List<String> args = getInputFiles();
			
			PrintWriter pw = null;
			SamReader sfr=null;
			try
				{
				this.file_index=0;
				
				if(args.isEmpty())
					{
					this.filenames.add("stdin");
					sfr=openSamReader(null);
					run(sfr);
					sfr.close();
					}
				else 
					{
					for(String filename: args)
						{
						filenames.add(filename);
						}
					
					for(String filename:filenames)
						{
						LOG.info("Opening "+filename);
						sfr=openSamReader(filename);
						run(sfr);
						sfr.close();
						this.file_index++;
						}
					}
				
				pw = openFileOrStdoutAsPrintWriter();
				
				if(this.print_header)
					{
					pw.print("#");
					if(this.print_chromstartend)
						{
						pw.print("chrom\tstart\tend\t");
						}
					pw.print("Transcript");
					for(String filename:filenames)
						{
						pw.print("\t"+filename);
						}
					pw.println();
					}
				
				List<Transcript> ordered=new ArrayList<Transcript>(this.name2transcript.values());
				Collections.sort(ordered,new java.util.Comparator<Transcript>()
						{
						@Override
						public int compare(Transcript a,Transcript b)
							{
							int i= a.tid-b.tid;
							if(i!=0) return i;
							i= a.start-b.start;
							if(i!=0) return i;
							i= a.end-b.end;
							if(i!=0) return i;
							return a.name.compareTo(b.name);
							}
						});
				
				for(Transcript tr:ordered)
					{
					if(tr.bad_flag) continue;
					
					if(this.removeZero)
						{
						//remove unmapped transcripts
						int x=0;
						for(x=0;x< tr.count.length;++x)
							if(tr.count[x]!=0) break;
						if(x==tr.count.length) continue;
						}
					if(this.print_chromstartend)
						{
						pw.print(firstDict.getSequence(tr.tid).getSequenceName());
						pw.print("\t");
						pw.print(tr.start);
						pw.print("\t");
						pw.print(tr.end);
						pw.print("\t");
						}
					pw.print(tr.name);
					for(int C:tr.count)
						{
						pw.print("\t");
						pw.print(C);
						}
					pw.println();
					}
			
			
				return RETURN_OK;
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(sfr);
				CloserUtil.close(pw);
				}
			}
		}
	public static void main(String[] args) {
		new HtSeqCount().instanceMainWithExit(args);
		}
	}
