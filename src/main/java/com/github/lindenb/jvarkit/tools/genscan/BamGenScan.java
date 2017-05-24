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
* July 2014: using a cigar string, remove bad quality cigars, secondaries

*/
package com.github.lindenb.jvarkit.tools.genscan;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.imageio.ImageIO;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;

/**
 * BamGenScan
 *
 */


@Program(name="bamgenscan",description="Paint a Genome Scan picture from a BAM")
public class BamGenScan extends AbstractGeneScan
	{
	private static final Logger LOG = Logger.build(BamGenScan.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	private List<Input> inputs=new ArrayList<Input>();
	
	/** a sample name and a filename */
	private static class Input
		{
		Sample sample;
		/** filename can contain more than one filename, can be colon-separated */
		String filenames;
		}
	
	@Override
	protected void drawPoints(Graphics2D g)
		{
		//loop over each input
		for(int input_id=0;input_id<this.inputs.size();++input_id)
			{		
			Input input=this.inputs.get(input_id);
			/** all readers for this sample */
			List<SamReader> readers=new ArrayList<SamReader>();
			for(String filename:input.filenames.split("[\\:]+"))
				{
				if(filename.isEmpty()) continue;
				LOG.info("Reading header for "+filename);
				SamReader sfr=SamFileReaderFactory.mewInstance().open(new File(filename));
				readers.add(sfr);
				
				}
			/* loop over each chromosome */
			for(ChromInfo ci:this.chromInfos)
				{
				if(!ci.visible) continue;
				LOG.info("["+input_id+"]alloc "+ci.dictSequenceLength+" for "+ci.getSequenceName());
				int count[]=new int[ci.dictSequenceLength];
				Arrays.fill(count, 0);
				
				
				Shape clip=g.getClip();
				g.setClip(new Rectangle2D.Double(
						ci.x,input.sample.y,
						ci.width,input.sample.height
						));
				
				for(SamReader sfr:readers) 
					{
					SAMRecordIterator sli=sfr.query(ci.getSequenceName(), 0, 0,true);
					
					//sli.setEmitUncoveredLoci(true);
					while(sli.hasNext())
						{
						SAMRecord rec=sli.next();
						if(rec.getReadUnmappedFlag()) continue;
						if(rec.getDuplicateReadFlag()) continue;
						if(rec.getReadFailsVendorQualityCheckFlag()) continue;
						if(rec.isSecondaryOrSupplementary()) continue;
						Cigar cigar=rec.getCigar();
						if(cigar==null) continue;
						int pos0=rec.getAlignmentStart() - 1;
						for(CigarElement ce:cigar.getCigarElements())
							{
							htsjdk.samtools.CigarOperator op=ce.getOperator();
							if(!op.consumesReferenceBases()) continue;
							if(op.consumesReadBases())
								{
								for(int L=0;L< ce.getLength();++L)
									{
									if(pos0<0 || pos0>=count.length) continue;
									count[pos0]++;
									++pos0;
									}
								}
							else
								{
								pos0+=ce.getLength();
								}
							}						
						}
					sli.close();
					}
				int win_size=(int)Math.ceil(ci.dictSequenceLength.doubleValue()/ci.width);
				if(win_size<0) win_size=1;
				Point2D.Double prev_points[]=new Point2D.Double[4];
				Point2D.Double curr_points[]=new Point2D.Double[4];
				for(int n=0;n+win_size <= count.length;n+=win_size)
					{
					MinMaxDouble minmax_val=new MinMaxDouble();
					double count_2=0;
					double total_2=0;
					double array_median[]=new double[win_size];
					for(int i=0;i< win_size;++i)
						{
						minmax_val.visit(count[i+n]);
						count_2 += count[i+n];
						total_2 += 1;
						array_median[i]=count[i+n];
						}
					if(total_2==0 || !minmax_val.isValid()) continue;
					Arrays.sort(array_median);
					
					for(int k=0;k< curr_points.length;++k)
						{
						curr_points[k]=new Point2D.Double();
						curr_points[k].x= ci.x+(((double)n/ci.minmaxBase.getMax())*ci.width);
						}
					
					curr_points[0].y=input.sample.y +input.sample.height- input.sample.height*this.minMaxY.getFraction((count_2)/total_2);
					curr_points[1].y=input.sample.y +input.sample.height- input.sample.height*this.minMaxY.getFraction(minmax_val.getMin());
					curr_points[2].y=input.sample.y +input.sample.height- input.sample.height*this.minMaxY.getFraction(minmax_val.getMax());
					curr_points[3].y=input.sample.y +input.sample.height- input.sample.height*this.minMaxY.getFraction(array_median[array_median.length/2]);
					
					
					for(int k=0;k< curr_points.length;++k)
						{
						if(prev_points[k]!=null)
							{
							switch(k)
								{
								case 0: g.setColor(Color.GREEN);break;
								case 1: g.setColor(Color.YELLOW);break;
								case 2: g.setColor(Color.RED);break;
								case 3: g.setColor(Color.BLUE);break;
								}
							g.draw(new Line2D.Double(prev_points[k],curr_points[k]));
							}
						prev_points[k]=curr_points[k];
						}
					}	
				g.setClip(clip);
				} //end for ci
			for(SamReader sfr:readers) CloserUtil.close(sfr);
			}//end for input
		}
	
	
	
	@Override
	protected List<ChromInfo> getChromInfos() {
		return this.chromInfos;
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			
			if(args.isEmpty())
				{
				LOG.error("illegal number of arguments");
				return -1;
				}
			SAMSequenceDictionary dict=null;
			for(int i=0;i<args.size();++i)
				{
				Input input=new Input();
				input.filenames=args.get(i);
				input.sample=new Sample();
				input.sample.name=input.filenames;
				input.sample.sample_id=super.samples.size();
				input.sample.minmax.setMin(0);
				input.sample.minmax.setMax(50);
				super.samples.add(input.sample);
				this.inputs.add(input);
				
				for(String filename:input.filenames.split("[\\:]+"))
					{
					if(filename.isEmpty()) continue;
					LOG.info("["+inputs.size()+"]Reading header for "+filename);
					SamReader sfr=SamFileReaderFactory.mewInstance().open(new File(filename));
					SAMFileHeader h=sfr.getFileHeader();
					sfr.close();
					SAMSequenceDictionary d=h.getSequenceDictionary();
					if(d==null)
						{
						LOG.error("Cannot get sequence dictionary for "+filename);
						return -1;
						}
					if(dict==null)
						{
						dict=d;
						}
					else if(!SequenceUtil.areSequenceDictionariesEqual(d, dict))
						{
						LOG.error("Sequence dictionaries are not the same.");
						return -1;
						}
					}
				
				}
			LOG.info("number of input:"+inputs.size());
			if(dict==null )
				{
				LOG.error("No dictionary found");
				return -1;
				}
		
			for(SAMSequenceRecord rec: dict.getSequences())
				{
				ChromInfo ci=new ChromInfo();
				ci.dictSequenceLength=rec.getSequenceLength();
				ci.sequenceName=rec.getSequenceName();
				ci.minmaxBase.setMin(0);
				ci.minmaxBase.setMax(ci.dictSequenceLength.doubleValue());
				
				ci.tid=chromInfos.size();
				if(ci.tid!=rec.getSequenceIndex()) throw new IllegalStateException("boum");
				chromInfos.add(ci);
				}
				
			
			
			BufferedImage img=makeImage();

			
			if(outputFile==null)
				{
				showGui(img);
				}
			else
				{
				ImageIO.write(img, "JPG", outputFile);
				}
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		
		}
	public static void main(String[] args)
		{
		new BamGenScan().instanceMainWithExit(args);
		}
	
	}
