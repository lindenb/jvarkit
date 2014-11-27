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
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.FilteringIterator;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SequenceUtil;

import java.awt.AlphaComposite;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import javax.imageio.ImageIO;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Hershey;
import com.github.lindenb.jvarkit.util.picard.MergingSamRecordIterator;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class BamCmpCoverage extends AbstractCommandLineProgram
	{
	private Map<String,Integer> sample2column=new HashMap<>();
	private int imgageSize=1000;
	private BufferedImage image=null;
	private double marginWidth=0;
	private double sampleWidth=0;
	private int minDepth=0;
	private int maxDepth=1000;
	
	
	private class Depth
		{
		int tid=0;
		int pos=0;
		int depths[]=new int[sample2column.size()];
		}
	
	private final Line2D.Double segment=new Line2D.Double();
	private void paint(Graphics2D g,final Depth depth)
		{
		
		for(int i=0;i< depth.depths.length;++i)
			{	
			
			if(depth.depths[i]< this.minDepth) continue;
			if(depth.depths[i]> this.maxDepth) continue;
			
			double x=this.marginWidth+i*sampleWidth;
			x+= this.sampleWidth*(depth.depths[i]-(double)this.minDepth)/((double)this.maxDepth-(double)this.minDepth);
			
			for(int j=0;j< depth.depths.length;++j)
				{	
				if(depth.depths[j]< this.minDepth) continue;
				if(depth.depths[j]> this.maxDepth) continue;
				double y=this.marginWidth+j*sampleWidth;
				y+=  this.sampleWidth -  this.sampleWidth*(depth.depths[j]-(double)this.minDepth)/((double)this.maxDepth-(double)this.minDepth);
				BamCmpCoverage.this.segment.x1= x;
				BamCmpCoverage.this.segment.y1= y;
				BamCmpCoverage.this.segment.x1= x+0.5;
				BamCmpCoverage.this.segment.y1= y-0.5;
				g.draw(BamCmpCoverage.this.segment);
				
				}
			}
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/BamCmpCoverage";
		}
	
	@Override
	public String getProgramDescription() {
		return "Creates an image for a comparative view of the depths sample vs sample.";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -o (file.png) output image name . Required");
		out.println(" -w (int) image width . Default: "+this.imgageSize);
		out.println(" -m (int) min depth . Default: "+this.minDepth);
		out.println(" -M (int) max depth . Default: "+this.maxDepth);
		out.println(" -r (chrom:start-end) region . Optional.");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		String regionStr=null;
		File imgOut=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:w:m:M:r:"))!=-1)
			{
			switch(c)
				{
				case 'o': imgOut=new File(opt.getOptArg());break;
				case 'w': this.imgageSize=Integer.parseInt(opt.getOptArg());break;
				case 'm': this.minDepth=Integer.parseInt(opt.getOptArg());break;
				case 'M': this.maxDepth=Integer.parseInt(opt.getOptArg());break;
				case 'r': regionStr=opt.getOptArg();break;
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
		if(imgOut==null)
			{
			error("output image file not defined");
			return -1;
			}
		
		if(this.imgageSize<1)
			{
			error("Bad image size:" +this.imgageSize);
			return -1;
			}
		
		if(this.minDepth<0)
			{
			error("Bad min depth : "+this.minDepth);
			return -1;
			}
		if(this.minDepth>=this.maxDepth)
			{
			error("Bad min<max depth : "+this.minDepth+"<"+this.maxDepth);
			return -1;
			}
		
		final SamRecordFilter filter=new SamRecordFilter()
			{
			@Override
			public boolean filterOut(SAMRecord first, SAMRecord second) {
				return filterOut(first);
				}
			
			@Override
			public boolean filterOut(SAMRecord rec) {
				if(rec.getReadUnmappedFlag()) return true;
				if(rec.isSecondaryOrSupplementary()) return true;
				if(rec.getDuplicateReadFlag())return true;
				if(rec.getNotPrimaryAlignmentFlag()) return true;
				if(rec.getReadFailsVendorQualityCheckFlag())return true;
				return false;
				}
			};
		Set<File> files=new HashSet<File>();
		try
			{
			SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			srf.disable(SamReaderFactory.Option.EAGERLY_DECODE);
			srf.disable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS);
			srf.disable(SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS);
			
			for(int i=opt.getOptInd();i< args.length;++i)
				{
				File f=new File(args[i]);
				if(f.getName().endsWith(".list"))
					{
					info("Reading BAM list from "+f);
					BufferedReader in=IOUtils.openFileForBufferedReading(f);
					String line;
					while((line=in.readLine())!=null)
						{
						if(line.trim().isEmpty() || line.startsWith("#")) continue;
						files.add(new File(line));
						}
					in.close();
					}
				else
					{
					files.add(f);
					}
				}
			if(files.isEmpty())
				{
				error("No BAM defined");
				return -1;
				}
			
			Comparator<SAMRecord> comparator=new Comparator<SAMRecord>()
				{
				@Override
				public int compare(SAMRecord samRecord1, SAMRecord samRecord2)
					{
					final int refIndex1 = samRecord1.getReferenceIndex();
			        final int refIndex2 = samRecord2.getReferenceIndex();
			        if (refIndex1 == -1) {
			            return (refIndex2 == -1? 0: 1);
			        } else if (refIndex2 == -1) {
			            return -1;
			        }
			        final int cmp = refIndex1 - refIndex2;
			        if (cmp != 0)
			        	{
			            return cmp;
			        	}
			        return samRecord1.getAlignmentStart() - samRecord2.getAlignmentStart();
					}
				};
			List<SamReader> readers=new ArrayList<SamReader>(files.size());
			List<CloseableIterator<SAMRecord>> iterators=new ArrayList<CloseableIterator<SAMRecord>>(files.size());
			

			
			Set<String> samples=new TreeSet<String>();
			SAMSequenceDictionary dict=null;
			//scan samples names
			for(File bamFile:files)
				{
				SamReader r= srf.open(bamFile);
				readers.add(r);
				
				SAMFileHeader h=r.getFileHeader();
				if(h.getSortOrder()!=SortOrder.coordinate)
					{
					r.close();
					error("file "+bamFile+" not sorted on coordinate");
					return -1;
					}
				if(dict==null)
					{
					dict=h.getSequenceDictionary();
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(dict,h.getSequenceDictionary()))
					{
					error("Found more than one dictint sequence dictionary");
					return -1;
					}

				
				for(SAMReadGroupRecord rg:h.getReadGroups())
					{
					String sample=rg.getSample();
					if(sample==null) continue;
					samples.add(sample);
					}
				
				if(regionStr==null)
					{
					iterators.add(new FilteringIterator(r.iterator(),filter));
					}
				else
					{
					int colon=regionStr.indexOf(':');
					String chrom;
					int chromStart;
					int chromEnd;

					if(colon==-1)
						{
						chrom=regionStr;
						}
					else
						{
						chrom=regionStr.substring(0,colon);
						}
					
					SAMSequenceRecord ssr= dict.getSequence(chrom);
					if(ssr==null)
						{
						error("Chromosome "+chrom+" not present in dictionary");
						return -1;
						}
					int hyphen=regionStr.indexOf('-', colon+1);
					if(hyphen!=-1)
						{
						chromStart=Integer.parseInt(regionStr.substring(colon+1,hyphen));
						chromEnd=Integer.parseInt(regionStr.substring(hyphen+1));
						}
					else
						{
						chromStart = 0;
						chromEnd = ssr.getSequenceLength()-1;
						}
					if(chromStart<0 || chromEnd<chromStart)
						{
						error("bad position in "+regionStr);
						return -1;
						}
					iterators.add(new FilteringIterator(r.query(chrom,chromStart,chromEnd,false),filter));
					
					}
				}
			info("Samples:"+samples.size());
			for(String sample:samples)
				{
				this.sample2column.put(sample, this.sample2column.size());
				}
			
			//create merging sam-reader
			MergingSamRecordIterator iter=new MergingSamRecordIterator(comparator,iterators);

			
			//create image
			info("Creating image "+this.imgageSize+"x"+this.imgageSize);
			this.image=new BufferedImage(this.imgageSize, this.imgageSize, BufferedImage.TYPE_INT_RGB);
			Graphics2D g=this.image.createGraphics();
			this.marginWidth=this.imgageSize*0.05;
			double drawingWidth=(this.imgageSize-1)-marginWidth;
			this.sampleWidth=drawingWidth/samples.size();
			//g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, this.imgageSize, this.imgageSize);
			g.setColor(Color.BLACK);
			Hershey hershey =new Hershey();
			for(String sample_x:samples)
				{
				g.setColor(Color.BLACK);
				hershey.paint(g,
						sample_x,
						marginWidth + sample2column.get(sample_x)*sampleWidth,
						0.0,
						sampleWidth*0.9,
						marginWidth
						);
				
        		AffineTransform old=g.getTransform();
        		AffineTransform tr= AffineTransform.getTranslateInstance(
        				marginWidth,
        				marginWidth + sample2column.get(sample_x)*sampleWidth
        				);
        		tr.rotate(Math.PI/2);
        		g.setTransform(tr);
        		hershey.paint(g,
						sample_x,
						0.0,
						0.0,
						sampleWidth*0.9,
						marginWidth
						);        		//g.drawString(this.tabixFile.getFile().getName(),0,0);
        		g.setTransform(old);
				
				for(String sample_y:samples)
					{
					
					Rectangle2D rect=new Rectangle2D.Double(
							marginWidth + sample2column.get(sample_x)*sampleWidth,
							marginWidth + sample2column.get(sample_y)*sampleWidth,
							sampleWidth,
							sampleWidth
							);
					g.setColor(Color.GRAY);
					g.draw(new Line2D.Double(rect.getMinX(),rect.getMaxY(),rect.getMaxX(),rect.getMinY()));
					g.setColor(Color.BLACK);
					g.draw(rect);
					}
				}
			
			Composite oldComposite =g.getComposite();
			g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.3f));
			//preivous chrom
			int prev_tid=-1;
			ArrayList<Depth> depthList=new ArrayList<Depth>();
			Stroke oldStroke=g.getStroke();
			g.setStroke(new BasicStroke(0.5f));
			g.setColor(Color.BLACK);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
			info("Scanning bams...");
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				if(filter.filterOut(rec)) continue;
				progress.watch(rec);

				SAMReadGroupRecord gr=rec.getReadGroup();
				if(gr==null) continue;
				String sample=gr.getSample();
				if(sample==null) continue;
				int sample_id= this.sample2column.get(sample);
				
				Cigar cigar= rec.getCigar();
				if(cigar==null) continue;
				int refPos=rec.getAlignmentStart();
				while(!depthList.isEmpty())
					{
					Depth front=depthList.get(0);
					if(front.tid!=rec.getReferenceIndex() || front.pos < refPos )
						{
						paint(g,front);
						depthList.remove(0);
						continue;
						}
					else
						{
						break;
						}		
					}
				if(prev_tid!=rec.getReferenceIndex())
					{
					g.setColor(Color.BLACK);
					String chrom=rec.getReferenceName();
					if(chrom.toLowerCase().contains("x"))
						{
						g.setColor(Color.RED);
						}
					else if(chrom.toLowerCase().contains("y"))
						{
						g.setColor(Color.PINK);
						}
					prev_tid = rec.getReferenceIndex();
					}
				for(CigarElement ce:cigar.getCigarElements())
					{
					CigarOperator op=ce.getOperator();
					if(!op.consumesReferenceBases()) continue;
					if(op.consumesReadBases())
						{
						for(int i=0;i< ce.getLength();++i)
							{
							Depth depth=null;
							int pos=refPos+i;
							if(depthList.isEmpty())
								{
								depth=new Depth();
								depth.pos=pos;
								depth.tid=prev_tid;
								depthList.add(depth);
								}
							else if(depthList.get(depthList.size()-1).pos< pos)
								{
								Depth prev=depthList.get(depthList.size()-1);
								while(prev.pos< pos)
									{
									depth=new Depth();
									depth.pos=prev.pos+1;
									depth.tid=prev_tid;
									depthList.add(depth);
									prev=depth;
									}
								depth=prev;
								}
							else
								{
								int lastPos=depthList.get(depthList.size()-1).pos;
								int distance= lastPos-pos;

								depth = depthList.get((depthList.size()-1)-(distance));
								if(depth.pos!=pos)
									{
									error(" "+pos+" vs "+depth.pos+" "+lastPos);
									System.exit(-1);
									}
								
								}
							depth.depths[sample_id]++;
							}
						}
					refPos+=ce.getLength();
					}
				}
			while(!depthList.isEmpty())
				{
				paint(g,depthList.remove(0));
				}
			progress.finish();
			iter.close();
			g.setStroke(oldStroke);
			g.setComposite(oldComposite);
			g.dispose();
			//close readers
			for(SamReader r:readers) r.close();
			
			//save file
			if(imgOut.getName().toLowerCase().endsWith(".png"))
				{
				ImageIO.write(this.image, "PNG", imgOut);
				}
			else
				{
				ImageIO.write(this.image, "JPG", imgOut);
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
			
			}
		}
	public static void main(String[] args) {
		new BamCmpCoverage().instanceMainWithExit(args);
	}
	}
