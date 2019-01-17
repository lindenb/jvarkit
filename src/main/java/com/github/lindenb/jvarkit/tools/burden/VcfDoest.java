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


History:

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.ucsc.TabixKnownGeneFileReader;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**

BEGIN_DOC


END_DOC
*/

@Program(name="vcfdoest",
	description="generate Transcript information for DOEST test",
	keywords={"vcf","burden","doest"}
		)
public class VcfDoest
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfDoest.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-k","--kg"},description=KnownGene.OPT_KNOWNGENE_DESC+". File must be compressed and indexed with tabix.",required=true)
	private String knownGeneURI = KnownGene.getDefaultUri();

	@Parameter(names={"-knc","--keepnoncoding"},description="keep non coding transcripts")
	private boolean keepNonCoding = false;

	@Parameter(names={"-f","--function"},description="User defined R function to be called after each VCF")
	private String userDefinedFunName = "";
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	private TabixKnownGeneFileReader knownGenesTabix=null;
	
	private static class TranscriptInfo implements Comparable<TranscriptInfo>
		{
		String contig;
		int txStart;
		int txEnd;
		String transcriptName;
		byte strand;
		int exonCount;
		int transcriptLength;
		String archetypeName;
		int archetypeLength;
		int ctxStart;
		Allele ref;
		Allele alt;
		double maf=MafCalculator.NODATA;
		int indexInTranscript0;
		String overlapName;
		byte genotypes[];

		
		
		@Override
		public int compareTo(final TranscriptInfo o)
			{
			int i= contig.compareTo(o.contig);
			if(i!=0) return i;
			i= transcriptName.compareTo(o.transcriptName);
			if(i!=0) return i;
			i= ctxStart - o.ctxStart;
			if(i!=0) return i;
			i= ref.compareTo(o.ref);
			if(i!=0) return i;
			i= alt.compareTo(o.alt);
			if(i!=0) return i;
			return 0;
			}
		}
	
	private static class TranscriptInfoCmp implements Comparator<TranscriptInfo> {
		@Override
		public int compare(TranscriptInfo o1, TranscriptInfo o2) {
			return o1.compareTo(o2);
			}
		}
		
	private static class TranscriptInfoCodec extends AbstractDataCodec<TranscriptInfo> {
		@Override
		public TranscriptInfo decode(DataInputStream dis) throws IOException {
			final TranscriptInfo o=new TranscriptInfo();
			try {
				o.contig = dis.readUTF();
			} catch (Exception e) {
				return null;
				}
			o.txStart = dis.readInt();
			o.txEnd = dis.readInt();
			o.transcriptName = dis.readUTF();
			o.strand = dis.readByte();
			o.exonCount = dis.readInt();
			o.transcriptLength = dis.readInt();
			o.archetypeName = dis.readUTF();
			o.archetypeLength = dis.readInt();
			
			o.ctxStart = dis.readInt();
			o.ref= Allele.create(TranscriptInfoCodec.readString(dis), true);
			o.alt= Allele.create(TranscriptInfoCodec.readString(dis), false);
			o.maf = dis.readDouble();
			
			o.indexInTranscript0 =  dis.readInt();
			o.overlapName = dis.readUTF();
			
			int n = dis.readInt();
		
			o.genotypes = new byte[n];
			for(int i=0;i<n;++i) 
				{
				o.genotypes[i] = dis.readByte();
				}

			return o;
			}
		
		@Override
		public void encode(final DataOutputStream dos, final TranscriptInfo o) throws IOException {
			dos.writeUTF(o.contig);
			dos.writeInt(o.txStart);
			dos.writeInt(o.txEnd);
			dos.writeUTF(o.transcriptName);
			dos.writeByte(o.strand);
			dos.writeInt(o.exonCount);
			dos.writeInt(o.transcriptLength);
			dos.writeInt(o.txStart);

			dos.writeUTF(o.archetypeName);
			dos.writeInt(o.archetypeLength);
			
			dos.writeInt(o.ctxStart);
			TranscriptInfoCodec.writeString(dos, o.ref.getDisplayString());
			TranscriptInfoCodec.writeString(dos, o.alt.getDisplayString());
			dos.writeDouble(o.maf);
			dos.writeInt(o.indexInTranscript0);
			
			dos.writeUTF(o.overlapName);
			
			dos.writeInt(o.genotypes.length);
			for(byte g : o.genotypes) {
				dos.writeByte(g);
			}
		
			}
		
		@Override
		public AbstractDataCodec<TranscriptInfo> clone() {
			return new TranscriptInfoCodec();
			}
		}
	
	
	public VcfDoest()
		{
		}
	 

	private List<KnownGene> overlap(final Interval interval) {
		final List<KnownGene> genes= new ArrayList<>();
		Iterator<KnownGene> iter= this.knownGenesTabix.iterator(interval);
		while(iter.hasNext())
			{
			final KnownGene  kg= iter.next();
			if(kg.isNonCoding() && !this.keepNonCoding) continue;
			genes.add(kg);
			}
		return genes;
		}
	
	

	
	private void run(
			final LineIterator lr,
			final PrintWriter pw
			) throws IOException
		
		{
		SortingCollection<TranscriptInfo> sorting = null;
		CloseableIterator<TranscriptInfo> iter2=null;
		try {
			while(lr.hasNext()) {
				VCFIterator in = VCFUtils.createVCFIteratorFromLineIterator(lr,true);
				final VCFHeader header=in.getHeader();
				
				final Pedigree pedigree = Pedigree.newParser().parse(header);
				if(pedigree.isEmpty())
					{
					throw new IOException("No pedigree found in header VCF header. use VcfInjectPedigree to add it");
					}
		
				final SortedSet<Pedigree.Person> individuals = new TreeSet<>();
		
				for(final Pedigree.Person individual:pedigree.getPersons()) {
					if(individual.isAffected() || individual.isUnaffected())
						{
						individuals.add(individual);
						}
					}
				boolean first=true;
				pw.println("# samples ( 0: unaffected 1:affected)");
				pw.print("population <- data.frame(family=c(");			
				first=true;for(final Pedigree.Person person : individuals) { if(!first) pw.print(","); pw.print("\""+person.getFamily().getId()+"\"");first=false;}
				pw.print("),name=c(");
				first=true;for(final Pedigree.Person person : individuals) { if(!first) pw.print(","); pw.print("\""+person.getId()+"\"");first=false;}
				pw.print("),status=c(");
				first=true;for(final Pedigree.Person person : individuals) { if(!first) pw.print(","); pw.print(person.isUnaffected()?0:1);first=false;}
				pw.println("))");
				
				sorting = SortingCollection.newInstance(TranscriptInfo.class,
						new TranscriptInfoCodec(),
						new TranscriptInfoCmp(),
						this.writingSortingCollection.getMaxRecordsInRam(),
						this.writingSortingCollection.getTmpPaths()
						);
				
				sorting.setDestructiveIteration(true);

		
				final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			
				/** loop over variants */
				while(in.hasNext() &&  !pw.checkError())
					{
					final VariantContext ctx = progess.watch(in.next());
					if(ctx.isFiltered()) continue;
					if(ctx.getAlternateAlleles().isEmpty()) continue;
					final Allele altAllele=ctx.getAltAlleleWithHighestAlleleCount();
					
					
					final MafCalculator mafCalculator= new MafCalculator(altAllele,ctx.getContig());
					boolean genotyped=false;
					for(final Pedigree.Person p:pedigree.getPersons())
						{
						if(!(p.isAffected() || p.isUnaffected())) continue;
						final Genotype g= ctx.getGenotype(p.getId());
						if(g==null) throw new IOException("Strange I cannot find individual "+p + " in the pedigree. Aborting.");
						if(g.isCalled()) {
							mafCalculator.add(g, p.isMale());
						}
						
						if(g.isHet() || g.isHomVar()) {
							if(!g.getAlleles().contains(altAllele)) continue;
							genotyped=true;
							break;
							}
						}
					if(!genotyped) continue;
					
					
					
					
					final Interval interval= new Interval(
							ctx.getContig(), 
							ctx.getStart(),
							ctx.getEnd()
							);
					
					final List<KnownGene> genes = this.overlap(interval);
					if(genes.isEmpty())  continue;
					/* find gene archetype = longest */
					
					for(final KnownGene kg:genes) {
						final TranscriptInfo trInfo = new TranscriptInfo();
						trInfo.contig = kg.getContig();
						trInfo.txStart = kg.getTxStart();
						trInfo.txEnd = kg.getTxEnd();
						trInfo.transcriptName = kg.getName();
						trInfo.strand = (byte)(kg.isPositiveStrand()?'+':'-');
						trInfo.exonCount = kg.getExonCount();
						trInfo.transcriptLength = kg.getTranscriptLength();
						
						trInfo.ctxStart = ctx.getStart();
						trInfo.ref = ctx.getReference();
						trInfo.alt = altAllele;
						trInfo.maf = mafCalculator.getMaf();
						
						trInfo.genotypes = new byte[individuals.size()];
						
						
						int idx=0;
						for(final Pedigree.Person individual:individuals) {
							final Genotype genotype=ctx.getGenotype(individual.getId());
							final byte b ;
							if (genotype.isHomRef()) {
								b=0;
							} else if (genotype.isHomVar() && genotype.getAlleles().contains(altAllele)) {
								b=2;
							} else if (genotype.isHet() && 
									genotype.getAlleles().contains(altAllele) &&
									genotype.getAlleles().contains(ctx.getReference())) {
								b=1;
							}
							/* we treat 0/2 has hom-ref */
							else if (genotype.isHet() && 
									!genotype.getAlleles().contains(altAllele) &&
									genotype.getAlleles().contains(ctx.getReference())) {
								LOG.warn("Treating "+genotype+" as hom-ref (0) alt="+altAllele);
								b=0;
								
							} 
							/* we treat 2/2 has hom-ref */
							else if (genotype.isHomVar() && !genotype.getAlleles().contains(altAllele)) {
								LOG.warn("Treating "+genotype+" as hom-ref (0) alt="+altAllele);
								b=0;
							}
							else {
								b=-9;
							}
							trInfo.genotypes[idx]=b;
							++idx;
						}
		
						
						
						
						KnownGene archetype=kg;
						/* find gene archetype = longest overlapping */
						for(final KnownGene kg2:genes) {
							if(kg2==kg) continue;
							if(archetype.getStrand().equals(kg2.getStrand()) && archetype.getTranscriptLength()< kg2.getTranscriptLength()) {
								archetype=kg2;
							}
						}
						
						trInfo.archetypeName = archetype.getName();
						trInfo.archetypeLength = archetype.getTranscriptLength();
						
						
						
						boolean ctxWasFoundInExon=false;
						final int ctxPos0=ctx.getStart()-1;
						int indexInTranscript0=0;
						
						for(final KnownGene.Exon exon:kg.getExons()) {
							//variant in exon ?
							if(!(exon.getStart()>(ctx.getEnd()-1)  || (ctx.getStart()-1)>=exon.getEnd())) {
								ctxWasFoundInExon=true;
								indexInTranscript0+=(ctxPos0-exon.getStart());
								
								if(kg.isNegativeStrand()) {
									indexInTranscript0 =  (kg.getTranscriptLength()-1)-indexInTranscript0;
								}
								trInfo.indexInTranscript0 = indexInTranscript0;
								trInfo.overlapName = exon.getName();
								sorting.add(trInfo);
								break;
							} else
							{
								indexInTranscript0+= (exon.getEnd()-exon.getStart());
							}
						}
					
					if(ctxWasFoundInExon) {
						continue;
					}
					
					
					indexInTranscript0=0;
					//search closest intron/exon junction
					for(int ex=0;ex+1<kg.getExonCount();++ex) {
						final KnownGene.Exon exon1= kg.getExon(ex  );
						indexInTranscript0+= (exon1.getEnd()-exon1.getStart());
						
						final KnownGene.Exon exon2= kg.getExon(ex+1);
						if(exon1.getEnd()<=ctxPos0 && ctxPos0<exon2.getStart())
							{
							final int dist_to_exon1 = ctxPos0 - exon1.getEnd();
							final int dist_to_exon2 = exon2.getStart() - ctxPos0;
							if(dist_to_exon2<dist_to_exon1) {
								indexInTranscript0++;
							}
							
							
							if(kg.isNegativeStrand()) {
								indexInTranscript0 =  (kg.getTranscriptLength()-1)-indexInTranscript0;
								}
							trInfo.indexInTranscript0=indexInTranscript0;
							trInfo.overlapName = exon1.getNextIntron().getName();
							sorting.add(trInfo);
							break;
							}
						}
					} // end loop over genes
				}//end while loop over variants
			progess.finish();
			sorting.doneAdding();
			LOG.info("done adding");
			iter2 = sorting.iterator();
			final EqualRangeIterator<TranscriptInfo> eqiter = new EqualRangeIterator<TranscriptInfo>(
					iter2,
					new Comparator<TranscriptInfo>() {
						@Override
						public int compare(final TranscriptInfo o1,final TranscriptInfo o2) {
							int i= o1.contig.compareTo(o2.contig);
							if(i!=0) return i;
							i= o1.transcriptName.compareTo(o2.transcriptName);
							return i;
							}
						}
					);
			while(eqiter.hasNext()) {
				final List<TranscriptInfo> list = eqiter.next();
				final TranscriptInfo front = list.get(0);
				
				pw.println("# BEGIN TRANSCRIPT "+front.transcriptName+" ##########################################");
	
				
				pw.println("transcript.chrom <- \""+front.contig+"\"");
				pw.println("transcript.txStart0 <- "+front.txStart+"");
				pw.println("transcript.txEnd0 <- "+front.txEnd+"");
				pw.println("transcript.name <- \""+front.transcriptName+"\"");
				pw.println("transcript.strand <- \""+((char)front.strand)+"\"");
				pw.println("transcript.length <- "+front.transcriptLength+"");
				pw.println("transcript.exonCount <- "+front.exonCount+"");
				pw.println("archetype.name <- \""+front.archetypeName+"\"");
				pw.println("archetype.length <- "+front.archetypeLength+"");
	
				
				pw.print("variants <- data.frame(chrom=c(");
				first=true; for(final TranscriptInfo v: list) {if(!first) pw.print(",");pw.print("\""+v.contig+"\"");first=false;}
				pw.print("),chromStart=c(");
				first=true; for(final TranscriptInfo v: list) {if(!first) pw.print(",");pw.print(v.ctxStart);first=false;}
				pw.print("),chromEnd=c(");
				first=true; for(final TranscriptInfo v: list) {if(!first) pw.print(",");pw.print(v.ctxStart+v.ref.length()-1);first=false;}
				pw.print("),refAllele=c(");
				first=true; for(final TranscriptInfo v: list) {if(!first) pw.print(",");pw.print("\""+v.ref.getDisplayString()+"\"");first=false;}
				pw.print("),altAllele=c(");
				first=true; for(final TranscriptInfo v: list) {if(!first) pw.print(",");pw.print("\""+v.alt.getDisplayString()+"\"");first=false;}
				pw.print("),positionInTranscript1=c(");
				first=true; for(final TranscriptInfo v: list) {if(!first) pw.print(",");pw.print(v.indexInTranscript0+1);first=false;}
				pw.print("),maf=c(");
				first=true; for(final TranscriptInfo v: list) {if(!first) pw.print(",");pw.print(v.maf);first=false;}
				pw.print("),overlapName=c(");
				first=true; for(final TranscriptInfo v: list) {if(!first) pw.print(",");pw.print("\""+v.overlapName+"\"");first=false;}
				pw.println("))");
	
				pw.println("# genotypes as a list. Should be a multiple of length(samples).");
				pw.println("# 0 is homref (0/0), 1 is het (0/1), 2 is homvar (1/1)");
				pw.println("# if the variant contains another ALT allele: (0/2) and (2/2) are considered 0 (homref)");
				pw.print("genotypes <- c(");
				first=true;
				for(final TranscriptInfo tr:list) {
					for(byte g: tr.genotypes) {
						if(!first) pw.print(","); first=false;
						pw.print((int)g);
					}
				}
				pw.println(")");
				pw.println("stopifnot(NROW(variants) * NROW(population) == length(genotypes) )");
	
				if(this.userDefinedFunName==null || this.userDefinedFunName.trim().isEmpty()) {
					pw.println("## WARNING not user-defined R function was defined");
					}
				else
					{
					pw.println("# consumme data with user-defined R function ");
					pw.println(this.userDefinedFunName+"()");
					}
				
				pw.println("# END TRANSCRIPT "+front.transcriptName+" ##########################################");
	
			} // end while eqiter
				
		
			eqiter.close();
			iter2.close();iter2=null;
				
			sorting.cleanup();sorting=null;
			}
		
		}
	finally 
		{
		CloserUtil.close(iter2);
		if(sorting!=null) sorting.cleanup();
		}
			
	}
	
	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.knownGeneURI==null || this.knownGeneURI.trim().isEmpty())
			{
			LOG.error("undefined option knownGeneURI");
			return -1;
			}
		if(!this.knownGeneURI.endsWith(".gz"))
			{
			LOG.error(this.knownGeneURI+" doesn't end with '.gz'");
			return -1;
			}
		LOG.info("reading "+this.knownGeneURI);
		BufferedReader r=null;
		LineIterator li=null;
		PrintWriter pw=null;
		try {
			LOG.info("loading tabix knownGene :"+this.knownGeneURI);
			this.knownGenesTabix = new TabixKnownGeneFileReader(this.knownGeneURI);
			final String inputFile = oneFileOrNull(args);
			
			pw=super.openFileOrStdoutAsPrintWriter(this.outputFile);
			
			li = (inputFile==null?
					IOUtils.openStreamForLineIterator(stdin()):
					IOUtils.openFileForLineIterator(new File(inputFile))
					);
			run(li,pw);
			CloserUtil.close(li);li=null;
			pw.flush();pw.close();pw=null;
			LOG.info("done");
			return RETURN_OK;
		} catch (Exception e) {
			LOG.error(e);
			return  -1;
		} finally 
		{
			CloserUtil.close(r);
			CloserUtil.close(li);
			CloserUtil.close(pw);
			CloserUtil.close(this.knownGenesTabix);
			
		}
		
		}
	 	
	
	public static void main(String[] args)
		{
		new VcfDoest().instanceMainWithExit(args);
		}
	}
