/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.bam2svg;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Hershey;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.ns.XLINK;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.svg.SVG;

public class BamToSVG extends AbstractCommandLineProgram
{
	private Hershey hershey=new Hershey();
	private IndexedFastaSequenceFile indexedFastaSequenceFile;
	private Interval interval;
	private Map<String, Sample> sampleHash=new HashMap<String, Sample>();
	private int drawinAreaWidth=1000;
	private boolean showClipping=false;
	private int featureHeight = 20;
	
	private static class Sample
		{
		String name;
		List<List<SAMRecord>> lines=new ArrayList<List<SAMRecord>>();
		}
	
	private static class Interval
		{
		String chrom;
		int start;
		int end;
		public String getChrom()
			{	
			return chrom;
			}
		public int getStart()
			{
			return start;
			}
		public int getEnd()
			{
			return end;
			}
		private int distance()
			{
			return this.getEnd()-this.getStart();
			}
		@Override
		public String toString() {
			return chrom+":"+start+"-"+end;
			}
		public boolean contains(int pos)
			{
			return start<=pos && pos<end;
			}
		
		}

	
	private BamToSVG()
		{
		}
	
	@Override
	public String getProgramDescription() {
		return "";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-R (file) "+getMessageBundle("reference.faidx"));
		out.println("-L (chrom:start-end) "+getMessageBundle("bed.interval0"));
		super.printOptions(out);
		}
	
	private int trim(int pos0)
		{
		return Math.min(Math.max(pos0, interval.start),interval.end);
		}

	
	private int left(final SAMRecord rec)
		{
		return (this.showClipping?rec.getUnclippedStart():rec.getAlignmentStart());
		}

	private int right(final SAMRecord rec)
		{
		return (this.showClipping?rec.getUnclippedEnd():rec.getAlignmentEnd());
		}
	private double baseToPixel(int pos)
		{
		return  ((pos - this.interval.getStart())/(double)this.interval.distance())*(this.drawinAreaWidth)
				;
		}

	
	private void readBamStream(SAMRecordIterator iter) throws IOException
		{
		while(iter.hasNext())
			{
			SAMRecord rec = iter.next();
			if(rec.getReadUnmappedFlag()) continue;
			if( !rec.getReferenceName().equals(this.interval.getChrom())) continue;
			if( right(rec)  < this.interval.getStart());
			if( left(rec)  >= this.interval.getEnd());
			String sampleName="";
			SAMReadGroupRecord srg = rec.getReadGroup();
			if(srg!=null)
				{
				String sm=srg.getSample();
				if(sm!=null) sampleName=sm;
				}
			Sample sample= this.sampleHash.get(sampleName);
			if(sample==null)
				{
				sample = new Sample();
				this.sampleHash.put(sampleName,sample);
				}
			
			int y=0;
			for(y=0;y< sample.lines.size();++y)
				{
				List<SAMRecord> line= sample.lines.get(y);
				final SAMRecord last = line.get(line.size()-1);
				if( right(last)+1 < left(rec) )
					{
					line.add(rec);
					break;
					}
				}
			if( y == sample.lines.size())
				{
				List<SAMRecord> line= new ArrayList<SAMRecord>();
				line.add(rec);
				sample.lines.add(line);
				}
			}
		}
	
	private void printGradientDef(
			XMLStreamWriter w,
			String gradId,
			String styleTop,//e.g: "stop-color:black;stop-opacity:1;"
			String styleMid //e.g: "stop-color:white;stop-opacity:1;"
			) throws XMLStreamException
		{
		w.writeStartElement("linearGradient");
		w.writeAttribute("id",gradId);
		w.writeAttribute("x1","50%");
		w.writeAttribute("x2","50%");
		w.writeAttribute("y1","0%");
		w.writeAttribute("y2","100%");
		w.writeEmptyElement("stop");
			w.writeAttribute("offset","0%");
			w.writeAttribute("style",styleTop);
		w.writeEmptyElement("stop");
			w.writeAttribute("offset","50%");
			w.writeAttribute("style",styleMid);
		w.writeEmptyElement("stop");
			w.writeAttribute("offset","100%");
			w.writeAttribute("style",styleTop);
		w.writeEndElement();
		}
	
	@Override
	public int doWork(String[] args)
		{
		GenomicSequence genomicSequence=null;
		String intervalStr=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"R:L:"))!=-1)
			{
			switch(c)
				{
				case 'L': intervalStr= opt.getOptArg();break;
				case 'R':
					{
					try
						{
						indexedFastaSequenceFile=new IndexedFastaSequenceFile(new File(opt.getOptArg())); 
						}
					catch(FileNotFoundException err)
						{
						error(err);
						return -1;
						}
					break;
					}
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		/* parse interval */
		if(intervalStr==null)
			{
			error(getMessageBundle("bed.interval0.undefined"));
			}
		int colon=intervalStr.indexOf(':');
		int hyphen=intervalStr.indexOf('-',colon+1);
		if(colon<1 || hyphen<=colon || hyphen+1==intervalStr.length())
			{
			System.err.println("Bad interval "+intervalStr);
			return -1;
			}
		
		this.interval=new Interval();
		this.interval.chrom=intervalStr.substring(0,colon);
		this.interval.start=Integer.parseInt(intervalStr.substring(colon+1,hyphen))+1;
		this.interval.end=Integer.parseInt(intervalStr.substring(hyphen+1));
			
		/* get genomic sequence */
		if(this.indexedFastaSequenceFile!=null)
			{
			genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, this.interval.chrom);
			}
		
		SamReader in=null;
		SAMRecordIterator iter=null;
		SamReaderFactory sfrf= SamReaderFactory.makeDefault();
		sfrf.validationStringency( ValidationStringency.SILENT);
		XMLStreamWriter w=null;
		try
			{
			/* read SAM data */
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				in = sfrf.open(SamInputResource.of(System.in));
				iter=in.iterator();
				readBamStream(iter);
				iter.close();
				in.close();
				}
			else
				{
				for(int i=opt.getOptInd();i< args.length;++i)
					{
					File filename=new File(args[i]);
					info("Reading from "+filename);
					in=sfrf.open(SamInputResource.of(filename));
					if(in.hasIndex())
						{
						iter=in.query(this.interval.getChrom(), this.interval.getStart(), this.interval.getEnd(), false);
						}
					else
						{
						info("Bam file not indexed !! "+filename);
						iter=in.iterator();
						}
					readBamStream(iter);
					iter.close();
					in.close();
					}
				}
			
			
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			w=xof.createXMLStreamWriter(System.out, "UTF-8");
			w.writeStartDocument("UTF-8", "1.0");
			w.writeStartElement("svg");
			w.writeDefaultNamespace(SVG.NS);
			w.writeNamespace("xlink", XLINK.NS);
			
			
			w.writeStartElement(SVG.NS,"title");
			w.writeCharacters(intervalStr);
			w.writeEndElement();
			
			w.writeStartElement(SVG.NS,"description");
			w.writeCharacters("Cmd:"+getProgramCommandLine()+"\n");
			w.writeCharacters("Version:"+getVersion()+"\n");
			w.writeCharacters("Author:"+getAuthorName()+" "+getAuthorMail()+"\n");
			w.writeCharacters("WWW:"+getOnlineDocUrl()+"\n");
			w.writeCharacters("Htsjdk: "+HtsjdkVersion.getHome()+" "+HtsjdkVersion.getVersion()+"\n");
			w.writeEndElement();
			
			
			w.writeStartElement(SVG.NS,"style");
			w.writeCharacters(
					".bA {}\n" + 
					".bC {}\n" +
					".bG {}\n" +
					".bT {}\n" +
					".bN {}\n" +
					""
					);
			w.writeEndElement();//style
			
			w.writeStartElement(SVG.NS,"defs");
			
			
			
			
			if(indexedFastaSequenceFile!=null)
				{
				for(String base:new String[]{"a","A","t","T","g","G","c","C","n","N"})
					{
					double width=drawinAreaWidth/(double)this.interval.distance();
					w.writeEmptyElement("path");
					w.writeAttribute("id","base"+base);
					w.writeAttribute("title",base);
					w.writeAttribute("class","b"+base.toUpperCase());
					w.writeAttribute("d",this.hershey.svgPath(
							base,
							0,
							0,
							width*0.95,
							featureHeight*0.95
							));
					}
				}
			
			
			w.writeEndElement();//defs
			int y=0;
			//loop over each sample
			for(Sample sample:this.sampleHash.values())
				{
				w.writeStartElement(SVG.NS,"g");
				w.writeAttribute("transform", "translate(0,"+y+")");
				
				/* print that line */
				for(int nLine=0;nLine< sample.lines.size();++nLine)
					{
					List<SAMRecord> line= sample.lines.get(nLine);
					double y_top= featureHeight*nLine;
					double mid_y= y_top+this.featureHeight/2.0;
					//loop over records on that line
					for(SAMRecord record: line)
						{
						w.writeStartElement(SVG.NS,"g");
						String title=record.getReadName();
						w.writeAttribute("title",title);
						
						/* print that sam record */
						final int unclipped_start= record.getUnclippedStart();
						Cigar cigar = record.getCigar();
						if(cigar==null) continue;
						byte bases[]=record.getReadBases();
						if(bases==null) continue;
						
						
						
						
						int readPos=0;
						Map<Integer,String> pos2insertions=new HashMap<Integer,String>();
						List<CigarElement> cigarElements= cigar.getCigarElements();
						

						
						int refPos=unclipped_start;
						/* loop over cigar string */
						for(int cidx=0; cidx< cigarElements.size(); cidx++ )
							{
							CigarElement ce = cigarElements.get(cidx);
							CigarOperator op=ce.getOperator();
							switch(ce.getOperator())
								{
								case D:
								case N:
									{
									int c_start = trim(refPos);
									int c_end   = trim(refPos + ce.getLength());
									if(c_start<c_end)
										{
										w.writeEmptyElement("line");
										w.writeEmptyElement("indel");
										w.writeAttribute("title", op.name()+ce.getLength());
										w.writeAttribute("x1", String.valueOf(baseToPixel(c_start)));
										w.writeAttribute("x2", String.valueOf(baseToPixel(c_end)));
										w.writeAttribute("y1", String.valueOf(mid_y));
										w.writeAttribute("y2", String.valueOf(mid_y));
										}
									
									refPos += ce.getLength();
									break;
									}
								case I: 
									{
									StringBuilder sb=new StringBuilder();
									for(int i=0;i< ce.getLength();++i)
										{
										sb.append((char)bases[readPos++]);
										}
									pos2insertions.put(refPos, sb.toString());
									break;
									}
								case H:
									{
									if(!this.showClipping)
										{
										refPos+=ce.getLength();
										break;
										}
									//NO break;
									}
								case S:
									{
									if(!this.showClipping)
										{
										readPos+=ce.getLength();
										refPos+=ce.getLength();
										break;
										}
									//NO break;
									}
								case X:
								case EQ:
								case M:
									{
									int match_start = refPos;
									int match_end = refPos + ce.getLength();
									
									//print sam background
									StringBuilder sb=new StringBuilder();
									if(record.getReadNegativeStrandFlag())
										{
										sb.append(" M ").append(baseToPixel(match_start)).append(',').append(y_top);
										sb.append(" h ").append(baseToPixel(match_end)-baseToPixel(match_start));
										sb.append(" v ").append(featureHeight);
										sb.append(" h ").append(-(baseToPixel(match_end)-baseToPixel(match_start)));
										sb.append(" Z");
										}
									else 
										{
										sb.append(" M ").append(baseToPixel(match_start)).append(',').append(y_top);
										sb.append(" h ").append(baseToPixel(match_end)-baseToPixel(match_start));
										sb.append(" v ").append(featureHeight);
										sb.append(" h ").append(-(baseToPixel(match_end)-baseToPixel(match_start)));
										sb.append(" Z");
										}
									w.writeEmptyElement("path");
									w.writeAttribute("d",sb.toString().trim());
									w.writeAttribute("class","r"+record.getFlags());
									
									if(op.consumesReadBases())
										{
										for(int i=0;i< ce.getLength();++i)
											{
											char ca=(char)bases[readPos];
											char cb='N';
											if(genomicSequence!=null)
												{
												cb=genomicSequence.charAt(refPos+i-1);
												}
											if(this.interval.contains(refPos+i))
												{
												w.writeEmptyElement("use");
												w.writeAttribute("x",String.valueOf( baseToPixel(refPos+i)));
												w.writeAttribute("y",String.valueOf(y_top));
												w.writeAttribute("ref", "base"+ca);
												}
											readPos++;
											}
										}
									

									
									refPos+=ce.getLength();
									break;
									}
								default:
									{
									throw new RuntimeException("Unknown SAM op: "+ce.getOperator());
									}
								}
							}	
						w.writeEndElement();//g
						}
					
					}
				
				/* surrounding frame for that sample */
				w.writeEmptyElement("rect");
				w.writeAttribute("class","frame");
				w.writeAttribute("x",String.valueOf(0));
				w.writeAttribute("y",String.valueOf(0));
				w.writeAttribute("width",String.valueOf(drawinAreaWidth));
				w.writeAttribute("height",String.valueOf(sample.lines.size()*featureHeight));//TODO

				
				y+= sample.lines.size()*featureHeight;//todo
				w.writeEndElement();//g for sample
				}
			
			
			w.writeEndElement();//svg
			w.writeEndDocument();
			w.flush();
			w.close();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(in);
			}
		}
	
	
	public static void main(String[] args)
		{
		new BamToSVG().instanceMainWithExit(args);
		}

}
