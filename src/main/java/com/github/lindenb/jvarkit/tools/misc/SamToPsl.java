/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.misc;

import java.awt.Color;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;


import htsjdk.samtools.SamReader;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.ucsc.PslAlign;

public class SamToPsl extends AbstractSamToPsl
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(AbstractSamToPsl.class);

	private PrintWriter out = null;
	private SamToPsl()
		{
		}
	
	
	
	private static class Align
		{
		int read0;
		int ref0;
		int len;
		Align(int read0,int ref0,int len)
			{
			this.read0=read0;
			this.ref0=ref0;
			this.len=len;
			}
		}
	private static List<Align> scancigar(int read0,int ref0,SAMRecord rec)
		{
		List<Align> L=new ArrayList<Align>();
		for(CigarElement ce:rec.getCigar().getCigarElements())
			{
			switch(ce.getOperator())
				{
				case I:
					{
					read0+=ce.getLength();
					break;
					}
				case N:
				case D:
					{
					ref0+=ce.getLength();
					break;
					}
				case M:
				case X:
				case EQ:
					{
					L.add(new Align(read0, ref0, ce.getLength()));
					read0+=ce.getLength();
					ref0+=ce.getLength();
					break;
					}
				case P:
				case H:
				case S:
					{
					break;
					}
				default:
					{
					throw new IllegalArgumentException("No handled:"+rec.getCigarString());
					}
				}
			}
		return L;
		}
	
	
	private List<PslAlign> makePslAlign(SAMRecord rec,SAMSequenceDictionary dict)
		{
		final List<PslAlign> aligns=new ArrayList<PslAlign>();

		final PslAlign a=new PslAlign();
		
		
		
		a.setStrand(rec.getReadNegativeStrandFlag()?'-':'+');
		a.setQName(rec.getReadName()+(rec.getReadPairedFlag()?(rec.getFirstOfPairFlag()?"/1":"/2"):"")+"_"+rec.getFlags());
		
		a.setTName(rec.getReferenceName());
		final SAMSequenceRecord ssr=dict.getSequence(rec.getReferenceIndex());
		if(ssr==null)
			{
			throw new RuntimeException("Cannot get SAMSequenceRecord dict for "+rec.getReferenceName());
			}
		a.setTSize( ssr.getSequenceLength());
		
		int qBaseInsert=0;
		int tBaseInsert=0;
		int readLength=0;
		for(final  CigarElement ce:rec.getCigar().getCigarElements())
			{
			switch(ce.getOperator())
				{
				case I: qBaseInsert+=ce.getLength();readLength+=ce.getLength();break;
				case N:
				case D: tBaseInsert+=ce.getLength();break;
				case M://threw
				case X:
				case EQ:
				case H:
				case S:
					{
					readLength+=ce.getLength();
					break;
					}
				default:break;
				}
			}
		a.setQBaseInsert(qBaseInsert);
		a.setTBaseInsert(tBaseInsert);
		
		boolean treat_as_paired=this.handle_paired_reads;
		if(rec.getSupplementaryAlignmentFlag()) treat_as_paired=false;

		if(rec.getReadPairedFlag())
			{
			if(rec.getMateUnmappedFlag())
				{
				treat_as_paired=false;
				}
			else
				{
				if(!rec.getMateReferenceIndex().equals(rec.getReferenceIndex()))
					{
					treat_as_paired=false;
					}
				else
					{
					if(rec.getReadNegativeStrandFlag()==rec.getMateNegativeStrandFlag())
						treat_as_paired=false;
					if(rec.getAlignmentStart() <= rec.getMateAlignmentStart())
						{
						if(rec.getMateAlignmentStart() <= rec.getAlignmentEnd())
							treat_as_paired=false;
						}
					}
				}
			}
		else
			{
			treat_as_paired=false;
			}
		
		
		if(!treat_as_paired)
			{
			a.setMatches(rec.getCigar().getReferenceLength());
			a.setQSize(readLength);
			a.setQEnd(readLength-(rec.getUnclippedEnd()-rec.getAlignmentEnd()));
			a.setTEnd(rec.getAlignmentEnd());
			
			final int readPos=rec.getAlignmentStart()-rec.getUnclippedStart();
			final int refPos=rec.getAlignmentStart()-1;
			a.setQStart(readPos);
			a.setTStart(refPos);
			
			for(final Align ctg:scancigar(readPos,refPos,rec))
				{
				a.addBlock(ctg.read0, ctg.ref0, ctg.len);
				}
			
			}
		else 
			{
			a.setMatches(rec.getCigar().getReferenceLength()+1);
			a.setQSize(readLength+1);
			
			
			
			if(rec.getAlignmentStart()< rec.getMateAlignmentStart())
				{
				
				a.setQEnd(readLength+1);
				a.setTEnd(rec.getMateAlignmentStart());
				
				int readPos=rec.getAlignmentStart()-rec.getUnclippedStart();
				int refPos=rec.getAlignmentStart()-1;
				a.setQStart(readPos);
				a.setTStart(refPos);
								
				for(final Align ctg:scancigar(readPos,refPos,rec))
					{
					a.addBlock(ctg.read0, ctg.ref0, ctg.len);
					}
				
				a.addBlock(
						readLength,
						rec.getMateAlignmentStart()-1,
						1
						);
				}
			else
				{
				a.setQStart(0);
				a.setQEnd((readLength+1)-(rec.getUnclippedEnd()-rec.getAlignmentEnd()));
				a.setTStart(rec.getMateAlignmentStart()-1);
				a.setTEnd(rec.getAlignmentEnd());
				
				
				a.addBlock(
						0,
						rec.getMateAlignmentStart()-1,
						1
						);
				
				for(final Align ctg:scancigar(1,rec.getAlignmentStart(),rec))
					{
					a.addBlock(ctg.read0, ctg.ref0, ctg.len);
					}
				}
			}
		
		aligns.add(a);
		
		return aligns;
		}
	
	private void scan(final SamReader in) 
		{
		final SAMSequenceDictionary dict=in.getFileHeader().getSequenceDictionary();
		if(dict==null) throw new RuntimeException("Sequence dictionary missing...");
		final SAMRecordIterator iter=in.iterator();
		final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);		
		while(iter.hasNext() && !this.out.checkError())
			{
			final  SAMRecord rec=progress.watch(iter.next());
			if(rec.getReadUnmappedFlag()) continue;
			for(final PslAlign a:makePslAlign(rec,dict))
				{
				out.println(toString(a,rec));
				}
			}
		progress.finish();
		iter.close();
		}
	
	private String toString(final PslAlign a,final  SAMRecord rec)
		{
		if(super.output_bed12)
			{
			final StringBuilder b=new StringBuilder();
		    //chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
			b.append(a.getTName());
			//chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
			b.append('\t');
			b.append(a.getTStart());
			//chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. 
			b.append('\t');
			b.append(a.getTEnd());
			//name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
			b.append('\t');
			b.append(a.getTName());
			//score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray: shade 	  	  	  	  	  	  	  	  	 
			b.append('\t');
			b.append((int)(((rec.getMappingQuality()>=255?0:rec.getMappingQuality())/255.0))*1000);
			//strand - Defines the strand - either '+' or '-'.
			b.append('\t');
			b.append((char)a.getStrand());
			//thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays).
			b.append('\t');
			b.append(a.getTStart());
			//thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
			b.append('\t');
			b.append(a.getTEnd());
			//itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
			b.append('\t');
			Color c=Color.BLACK;
			if(rec.getReadPairedFlag())
				{
				c=Color.ORANGE;
				if(rec.getProperPairFlag()) c=Color.GREEN;
				}
			b.append(c.getRed()).append(",").append(c.getGreen()).append(",").append(c.getBlue());
			//blockCount - The number of blocks (exons) in the BED line.
			b.append('\t');
			b.append(a.getBlockCount());
			//blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
			b.append('\t');
			for(PslAlign.Block block:a.getBlocks())
				{
				b.append(block.size());
				b.append(',');
				}
			//blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount. 
			b.append('\t');
			for(PslAlign.Block block:a.getBlocks())
				{
				b.append(block.getTStart()-a.getTStart());
				b.append(',');
				}
			return b.toString();
			}
		else
			{
			return a.toString();
			}
		}
	

	@Override
	protected Collection<Throwable> call(final String inputName) throws Exception {
		
		SamReader sfr=null;
		try
			{
			sfr = super.openSamReader(inputName);
			this.out = super.openFileOrStdoutAsPrintWriter();
			scan(sfr);
			this.out.flush();
			LOG.info("done");
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(sfr);
			CloserUtil.close(this.out);
			this.out=null;
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new SamToPsl().instanceMainWithExit(args);
		}
	}
