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
package com.github.lindenb.jvarkit.tools.sam2tsv;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.tools.misc.IlluminaReadName;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.ProgressLoggerInterface;

@Program(name="prettysam",
description="Pretty SAM alignments",
keywords={"sam","bam",}
)

public class PrettySam extends Launcher {
	private static final Logger LOG = Logger.build(PrettySam.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private File refFile = null;
	
	



	public static class PrettySAMWriter implements SAMFileWriter
		{
		private final NumberFormat fmt = new DecimalFormat("#,###");
		private final PrintWriter pw ;
		private SAMFileHeader header = null;
		private long nLine=0;
		private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
		private GenomicSequence genomicSequence=null;
		
		
		public PrettySAMWriter(PrintWriter pw) {
			this.pw = pw;
			
			
		}
		
		private static class Base
			{
			char readbase;
			char readqual;
			int  readpos;
			char refbase;
			int  refpos;
			CigarOperator cigaroperator;
			}

		
		public PrettySAMWriter setReferenceFile(final File f) throws IOException{
			if(f==null) return null;
			CloserUtil.close(this.indexedFastaSequenceFile);
			this.genomicSequence=null;
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(f);
			return this;
			}
		
		private char getReferenceAt(final String contig,int refpos) {
			if(this.indexedFastaSequenceFile==null) return 'N';
			
			if(this.genomicSequence==null || !this.genomicSequence.getChrom().equals(contig))
				{
				this.genomicSequence = new GenomicSequence(this.indexedFastaSequenceFile, contig);
				}
			if((refpos-1)>=this.genomicSequence.length())return 'N';
			return this.genomicSequence.charAt(refpos-1);
			}
		
		@Override
		public void setProgressLogger(ProgressLoggerInterface progress) {
			
			}
		@Override
		public SAMFileHeader getFileHeader() {
			return header;
		}	
		
		public void writeHeader(final SAMFileHeader header) {
			this.header = header;
			pw.flush();
			}
		
		private void label(int cols,final String s) { pw.printf("%"+cols+"s : ",s); }
		
		@Override
		public void addAlignment(final SAMRecord rec) {
			final int margin1=19;
			final int margin2=margin1+5;
			
			++this.nLine;
			if(this.nLine>1) pw.println();
			pw.println(">>>>> "+nLine);
			label(margin1,"Read-Name");pw.println(rec.getReadName());
			
			 new IlluminaReadName.Parser().apply(rec.getReadName()).
			 	ifPresent(ilmn->{
					label(margin2,"Instrument");pw.println(ilmn.getInstrument());
					label(margin2,"Run");pw.println(ilmn.getRunId());
					label(margin2,"FlowCell");pw.println(ilmn.getFlowCell());
					label(margin2,"Lane");pw.println(ilmn.getLane());
					label(margin2,"Tile");pw.println(ilmn.getTile());
					label(margin2,"X");pw.println(ilmn.getX());
					label(margin2,"Y");pw.println(ilmn.getY());
					});
			label(margin1,"Flag");pw.println(rec.getFlags());
			for(final SAMFlag flg:SAMFlag.values())
				{
				if(!flg.isSet(rec.getFlags())) continue;
				label(margin2,flg.getLabel());pw.println(flg.intValue());
				}
			if(!rec.getReadUnmappedFlag()) {
				
				label(margin1,"MAPQ");
				pw.print(rec.getMappingQuality());
				if(rec.getMappingQuality()==SAMRecord.UNKNOWN_MAPPING_QUALITY)
					{
					pw.print(" (unknown)");
					}
				pw.println();
				
				label(margin1,"Contig");pw.println(rec.getReferenceName()+"  (index:"+rec.getReferenceIndex()+")");
				if(rec.getUnclippedStart()!=rec.getStart())
					{
					label(margin1,"Unclipped-Start");pw.println(this.fmt.format(rec.getUnclippedStart()));
					}
				label(margin1,"Start");pw.println(this.fmt.format(rec.getAlignmentStart()));
				label(margin1,"End");pw.println(this.fmt.format(rec.getAlignmentEnd()));
				if(rec.getUnclippedEnd()!=rec.getEnd())
					{
					label(margin1,"Unclipped-End");pw.println(this.fmt.format(rec.getUnclippedEnd()));
					}
				label(margin1,"Strand");pw.println(rec.getReadNegativeStrandFlag()?"-":"+");
				}
			if(rec.getReadPairedFlag() && !rec.getMateUnmappedFlag())
				{
				if(rec.getInferredInsertSize()!=0)
					{
					label(margin1,"Insert-Size");pw.println(this.fmt.format(rec.getInferredInsertSize()));
					}
				label(margin1,"Mate-Contig");pw.println(rec.getMateReferenceName()+"  (index:"+rec.getMateReferenceIndex()+")");
				label(margin1,"Mate-Start");pw.println(this.fmt.format(rec.getMateAlignmentStart()));
				label(margin1,"Mate-Strand");pw.println(rec.getMateNegativeStrandFlag()?"-":"+");
				}
			final Cigar cigar=rec.getCigar();
			if(cigar!=null && !cigar.isEmpty())
				{
				label(margin1,"Cigar");
				pw.println(rec.getCigarString()+" (N="+cigar.numCigarElements()+")");
				}
			
			final List<Base> align = new ArrayList<>();
			final String bases = rec.getReadString();
			final String quals = rec.getBaseQualityString();

			if(rec.getReadUnmappedFlag() || cigar==null || cigar.isEmpty())
				{
				for(int i=0;i< bases.length();i++)
					{
					final Base b=new Base();
					b.readpos=i;
					b.refbase = ' ';
					b.readpos = -1;
					b.cigaroperator = CigarOperator.P;
					b.readbase = (bases!=null && i>=0 && i<bases.length()?bases.charAt(i):'*');
					b.readqual = (quals!=null && i>=0 && i<quals.length()?quals.charAt(i):'?');
					align.add(b);
					}
				}
			else
				{
				int refpos=rec.getUnclippedStart();
				int readpos=0;
				if(cigar.numCigarElements()>1 && cigar.getCigarElement(0).getOperator()==CigarOperator.H)
					{
					readpos = rec.getUnclippedStart()-rec.getAlignmentStart();//WILL be negative !
					}
				
				
				final Function<Integer, Character> pos2base= (i)->(bases!=null && i>=0 && i<bases.length()?bases.charAt(i):'*');
				final Function<Integer, Character> pos2qual= (i)->(quals!=null && i>=0 && i<quals.length()?quals.charAt(i):'*');
				final Function<Integer, Character> pos2ref= (i)->getReferenceAt(rec.getContig(),i);
				
				
				for(final CigarElement ce:cigar)
					{
					final CigarOperator op = ce.getOperator();
					switch(op)
						{
						case P: break;
						case D: case N:
							{
							final Base b=new Base();
							b.cigaroperator = op;
							b.refpos  = refpos;
							b.refbase = pos2ref.apply(refpos);
							b.readbase = '^';
							b.readqual = '^';
							b.readpos = -1;
							align.add(b);
							refpos+=ce.getLength();
							break;
							}
						case X:case EQ:case M:case S:case H:case I:
							{
							for(int i=0;i< ce.getLength();++i)
								{
								final Base b=new Base();
								b.refpos=refpos;
								b.cigaroperator = op;
								b.readbase = pos2base.apply(readpos);
								b.readqual = pos2qual.apply(readpos);
								b.readpos=readpos;
								
								readpos++;//ok with 'H', because negative from beginning *
								if(op.equals(CigarOperator.I))
									{
									b.refpos = -1;
									b.refbase  = '^';
									}
								else
									{
									b.refpos= refpos;
									b.refbase  = pos2ref.apply(refpos);
									refpos++;
									}
								
								align.add(b);
								}
							break;
							}
						}
					
					}
				
				}
			
			if(!align.isEmpty())
				{
				label(margin1,"Sequence");
				pw.println();
				int x=0;
				final int FASTA_LEN=60;
				while(x<align.size())
					{
					for(int side=0;side<5;++side)
						{
						switch(side)
							{
							case 0: label(margin2,"Read ("+ this.fmt.format(align.stream().
									skip(x).
									mapToInt(B->B.readpos).
									filter(X -> X>=0).
									findFirst().orElse(-1))+")");
									break;
							case 4: label(margin2,"Qual");break;
							case 2:
								if(cigar==null || cigar.isEmpty()) continue;
								label(margin2,"Ref ("+this.fmt.format(align.stream().
										skip(x).
										mapToInt(B->B.refpos).
										filter(X -> X>0).
										findFirst().orElse(-1))+")");
								break;
							case 3:
								if(cigar==null || cigar.isEmpty()) continue;
								label(margin2,"Op");break;
							case 1:
								if(cigar==null || cigar.isEmpty()) continue;
								label(margin2,"Mid");break;
							default:break;
							}
						for(int y=0;y<FASTA_LEN && x+y<align.size();++y)
							{
							final Base base = align.get(x+y);
							if(y!=0 && y%10==0) pw.print(" ");
							final char c;
							switch(side)
								{
								case 0: c =  base.readbase;break;
								case 4: c =  base.readqual;break;
								case 2: c =  base.refbase;break;
								case 3: c = base.cigaroperator.name().charAt(0);break;
								case 1:
									c= base.cigaroperator.isAlignment()?
											(base.refbase==base.readbase?'|':' ')
											: ' ';
									break;
								default : c='?'; break;
								}
							pw.print(c);
							}
						pw.println();
						
						}
					x+=FASTA_LEN;
					pw.println();
					}
				}
			
			
			pw.println("<<<<< "+nLine);
			pw.flush();
			}
		
		@Override
		public void close() {
			pw.close();
			CloserUtil.close(this.indexedFastaSequenceFile);
			this.genomicSequence=null;
			}
		}
	@Override
	public int doWork(final List<String> args) {
		SamReader r = null;
		PrettySAMWriter out = null;
		CloseableIterator<SAMRecord> iter = null;
		try 
			{
			r= super.openSamReader(oneFileOrNull(args));
			out = new PrettySAMWriter(super.openFileOrStdoutAsPrintWriter(this.outputFile));
			out.writeHeader(r.getFileHeader());
			out.setReferenceFile(this.refFile);
			iter = r.iterator();
			while(iter.hasNext())
				{
				out.addAlignment(iter.next());
				}
			out.close();out=null;
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(r);
			CloserUtil.close(out);
			}
		}
	
	public static void main(final String[] args) {
		new PrettySam().instanceMainWithExit(args);
	}

}
