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
package com.github.lindenb.jvarkit.tools.misc;

import java.awt.Color;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;


import htsjdk.samtools.SamReader;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.ucsc.PslAlign;

/**

BEGIN_DOC

## Motivation

Convert **SAM/BAM** to **PSL** http://genome.ucsc.edu/FAQ/FAQformat.html#format2 or **BED12** .

Properly-paired reads are extended to the mate's position.

What ? **bamtobed** http://bedtools.readthedocs.org/en/latest/content/tools/bamtobed.html does the same job ?! too late.

## Cited in:

   * "R2C2: Improving nanopore read accuracy enables the sequencing of highly-multiplexed full-length single-cell cDNA" biorxiv  https://doi.org/10.1101/338020 

### Example

```
$ samtools view -b  http://hgdownload-test.cse.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeq10t12C3hFR2x75Th131Il200AlnRep1.bam "chr15:81575506-81616397" |\
   java -jar ~/src/jvarkit-git/dist/sam2psl.jar -s   > out.psl

$ tail out.psl
6065	0	0	0	0	0	0	5964	-	HWI-ST0787:100:C02F9ACXX:3:1105:11756:193850_2:N:0:/2_147	101	0	101	chr15	10349497481616327	81622456	3	1,5,95,	0,1,6,	81616327,81616393,81622362,
6065	0	0	0	0	0	0	5964	-	HWI-ST0787:100:C02F9ACXX:3:1301:4643:94800_2:N:0:/2_147	101	0	101	chr15	103494974	81616334	81622456	3	1,5,95,	0,1,6,	81616334,81616393,81622362,
6065	0	0	0	0	0	0	5964	-	HWI-ST0787:100:C02F9ACXX:3:1308:18580:185579_2:Y:0:/2_147	101	0	101	chr15	10349497481616327	81622456	3	1,5,95,	0,1,6,	81616327,81616393,81622362,
6065	0	0	0	0	0	0	5964	-	HWI-ST0787:100:C02F9ACXX:3:2205:10117:76559_2:N:0:/2_147	101	0	101	chr15	10349497481616321	81622456	3	1,5,95,	0,1,6,	81616321,81616393,81622362,
6065	0	0	0	0	0	0	5964	-	HWI-ST0787:100:C02F9ACXX:3:2206:7885:15613_2:N:0:/2_403	101	0	101	chr15	103494974	81613633	81622456	3	1,5,95,	0,1,6,	81613633,81616393,81622362,
6065	0	0	0	0	0	0	5964	-	HWI-ST0787:100:C02F9ACXX:3:2206:7885:15613_2:N:0:/2_147	101	0	101	chr15	103494974	81616308	81622456	3	1,5,95,	0,1,6,	81616308,81616393,81622362,
6065	0	0	0	0	0	0	5964	-	HWI-ST0787:100:C02F9ACXX:3:2303:12879:149117_1:Y:0:/1_83	101	0	101	chr15	10349497481616334	81622456	3	1,5,95,	0,1,6,	81616334,81616393,81622362,
6064	0	0	0	0	0	0	5964	+	HWI-ST0787:100:C02F9ACXX:3:1301:11600:100190_1:N:0:/1_99	100	0	100	chr15	10349497481616393	81622457	2	4,96,	0,4,	81616393,81622361,
6064	0	0	0	0	0	0	5964	+	HWI-ST0787:100:C02F9ACXX:3:2304:5980:187674_1:Y:0:/1_99	100	0	100	chr15	103494974	81616393	81622457	2	4,96,	0,4,	81616393,81622361,
6065	0	0	0	0	0	0	5964	-	HWI-ST0787:100:C02F9ACXX:3:1306:18607:99733_2:N:0:/2_147	101	0	101	chr15	10349497481616334	81622457	3	1,4,96,	0,1,5,	81616334,81616394,81622362,
```

used as a custom track in the **UCSC genome browser**.

![img](http://i.imgur.com/Gi6Sd0M.png)


### See also

  * bedtools/bamtobed : http://bedtools.readthedocs.org/en/latest/content/tools/bamtobed.html

## Cited in

  *  "Depletion of hemoglobin transcripts and long read sequencing improves the transcriptome annotation of the polar bear (Ursus maritimus)
Ashley Byrne, Megan A Supple, Roger Volden, Kristin L Laidre, Beth Shapiro, Christopher Vollmers" bioRxiv 527978; doi: https://doi.org/10.1101/527978 

END_DOC
*/
@Program(name="sam2psl",
	deprecatedMsg="use bedtools/bamtobed",
	description="Convert SAM/BAM to PSL http://genome.ucsc.edu/FAQ/FAQformat.html#format2 or BED12",
	keywords={"sam","bam","psl"}
	)
public class SamToPsl extends Launcher
	{
	private static final Logger LOG = Logger.build(SamToPsl.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-s","--single"},description="treat all reads as single end.")
	private boolean handle_paired_reads=false;
	@Parameter(names={"-B","bed12"},description="Export as BED 12.")
	private boolean output_bed12=false;
	
	private PrintWriter out = null;
	
	
	
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
		final List<Align> L=new ArrayList<Align>();
		for(final CigarElement ce:rec.getCigar().getCigarElements())
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
			for(final PslAlign a: makePslAlign(rec,dict))
				{
				out.println(toString(a,rec));
				}
			}
		progress.finish();
		iter.close();
		}
	
	private String toString(final PslAlign a,final  SAMRecord rec)
		{
		if(this.output_bed12)
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
			b.append(a.getStrand());
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
	public int doWork(final List<String> args)
		{
		SamReader sfr=null;
		try
			{
			sfr = super.openSamReader(oneFileOrNull(args));
			this.out = super.openFileOrStdoutAsPrintWriter(outputFile);
			scan(sfr);
			this.out.flush();
			this.out.close();
			this.out = null;
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(sfr);
			CloserUtil.close(this.out);
			this.out=null;
			}
		}
	public static void main(final String[] args)
		{
		new SamToPsl().instanceMainWithExit(args);
		}
	}
