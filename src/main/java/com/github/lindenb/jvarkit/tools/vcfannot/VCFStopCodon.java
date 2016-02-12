/**
 * Author:
 * 	Pierre Lindenbaum PhD
 * Date:
 * 	Fev-2014
 * Contact:
 * 	plindenbaum@yahoo.fr
 * Motivation:
 * 	Idea from Solena: successive synonymous mutations are a stop codong
 */
package com.github.lindenb.jvarkit.tools.vcfannot;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;



/**
 * VCFStopCodon
 * @SolenaLS 's idea: consecutive synonymous bases give a stop codon
 *
 */
public class VCFStopCodon extends AbstractVCFStopCodon
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VCFStopCodon.class);

	private IntervalTreeMap<KnownGene> knownGenes=new IntervalTreeMap<KnownGene>();
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private GenomicSequence genomicSequence=null;
	private SortingCollection<Variant> variants= null;
	
	private static class MutedSequence extends DelegateCharSequence
		{
		private Map<Integer, Character> pos2char=new TreeMap<Integer, Character>();
		MutedSequence(CharSequence wild)
			{
			super(wild);
			}
		
		void put(int pos,char c)
			{
			if(pos2char.containsKey(pos)) throw new IllegalStateException();
			this.pos2char.put(pos, c);
			}
		
		@Override
		public char charAt(int i)
			{
			Character c= pos2char.get(i);
			return c==null?getDelegate().charAt(i):c;
			}
		
		}
	
	
	private static class ProteinCharSequence extends DelegateCharSequence
		{
		private GeneticCode geneticCode;
		ProteinCharSequence(GeneticCode geneticCode,CharSequence cDNA)
			{
			super(cDNA);
			this.geneticCode=geneticCode;
			}
		
		@Override
		public char charAt(int i)
			{
			return geneticCode.translate(
				getDelegate().charAt(i*3+0),
				getDelegate().charAt(i*3+1),
				getDelegate().charAt(i*3+2));
			}	
		
		@Override
		public int length()
			{
			return getDelegate().length()/3;
			}
		}
		
	private void loadKnownGenesFromUri() throws IOException
		{
		BufferedReader in = null;
		try {
			final SAMSequenceDictionary dict=this.indexedFastaSequenceFile.getSequenceDictionary();
	        if(dict==null) throw new IOException("dictionary missing");

			LOG.info("loading genes from "+this.kgURI);
			in =IOUtils.openURIForBufferedReading(this.kgURI);
			final Pattern tab=Pattern.compile("[\t]");
			String line = null;
			while((line=in.readLine())!=null)
				{
				line= line.trim();
				if(line.isEmpty()) continue;
				String tokens[]=tab.split(line);
				final KnownGene g=new KnownGene(tokens);
				if(g.isNonCoding()) continue;
				if(dict.getSequence(g.getContig())==null)
					{
					LOG.warn("Unknown chromosome "+g.getContig()+" in dictionary");
					continue;
					}
				//use 1 based interval
				final Interval interval=new Interval(g.getContig(), g.getTxStart()-1, g.getTxEnd());
				this.knownGenes.put(interval, g);
				}
			CloserUtil.close(in);in=null;
			LOG.info("genes:"+knownGenes.size());
			}
		finally
			{
			CloserUtil.close(in);
			}
		}
	
	
	
	private static int ID_GENERATOR=0;
	static private class Variant
		{
		String contig;
		int genomicPosition1=0;
		String transcriptName;
		Allele refAllele;
		Allele altAllele;
		int position_in_cdna=-1;
		String wildCodon=null;
		String mutCodon=null;
		int sorting_id;
		
		Variant()
			{
			
			}
		
		Variant(final VariantContext ctx,final Allele allele,final KnownGene gene) {
			this.contig = ctx.getContig();
			this.genomicPosition1=ctx.getStart();
			this.transcriptName = gene.getName();
			this.refAllele = ctx.getReference();
			this.altAllele = allele;
			}
		int positionInCodon() { return  position_in_cdna%3;}
		int codonStart() { return this.position_in_cdna - this.positionInCodon();}
		@Override
		public String toString() {
			return contig+"\t"+genomicPosition1+"\t"+refAllele.getBaseString()+"\t"+
					altAllele.getBaseString()+"\t"+
					transcriptName+"\t"+
					position_in_cdna+"\t"+
					codonStart()+"\t"+
					positionInCodon()+"\t"+
					wildCodon+"\t"+
					mutCodon
					;
			}
		
		}
	static private class VariantCodec extends AbstractDataCodec<Variant>
		{
		@Override
		public Variant decode(DataInputStream dis) throws IOException {
			String contig;
			try {
				contig = dis.readUTF();
			} catch (Exception e) {
				return null;
				}
			final Variant variant = new Variant();
			variant.contig = contig;
			variant.genomicPosition1 = dis.readInt();
			variant.transcriptName = dis.readUTF();
			variant.refAllele = Allele.create(dis.readUTF(), true);
			variant.altAllele = Allele.create(dis.readUTF(), false);
			variant.position_in_cdna = dis.readInt();
			variant.wildCodon = dis.readUTF();
			variant.mutCodon = dis.readUTF();
			variant.sorting_id = dis.readInt();
			return variant;
		}

		@Override
		public void encode(DataOutputStream dos, Variant v) throws IOException {
			dos.writeUTF(v.contig);
			dos.writeInt(v.genomicPosition1);
			dos.writeUTF(v.transcriptName);
			dos.writeUTF(v.refAllele.getBaseString());
			dos.writeUTF(v.altAllele.getBaseString());
			dos.writeInt(v.position_in_cdna);
			dos.writeUTF(v.wildCodon);
			dos.writeUTF(v.mutCodon);
			dos.writeInt(v.sorting_id);
		}

		@Override
		public VariantCodec clone() {
			return new VariantCodec();
		}
		
		}
	static private class VariantComparator implements Comparator<Variant>
		{
		SAMSequenceDictionary dict;
		VariantComparator(SAMSequenceDictionary dict) {
			this.dict = dict;
		}
		int contig(Variant v) { return dict.getSequenceIndex(v.contig);}
		@Override
		public int compare(Variant o1, Variant o2) {
			int i= contig(o1) - contig(o2);
			if(i!=0) return i;
			i= o1.transcriptName.compareTo(o2.transcriptName);
			if(i!=0) return i;
			i= o1.position_in_cdna-o2.position_in_cdna;
			if(i!=0) return i;
			return o1.sorting_id - o2.sorting_id;
			}
	
		}

	
	
	static final String TAG="STREAMCODON";
	
	
	
	protected Collection<Throwable> scan(final VcfIterator r) throws IOException {
			
			CloseableIterator<Variant> varIter = null;
			try {
			LOG.info("opening REF:"+referenceFile);
			
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(this.referenceFile);
	        final SAMSequenceDictionary dict=this.indexedFastaSequenceFile.getSequenceDictionary();
	        if(dict==null) throw new IOException("dictionary missing");
			loadKnownGenesFromUri();
			
			this.variants = SortingCollection.newInstance(Variant.class,
					new VariantCodec(),
					new VariantComparator(dict),
					super.getMaxRecordsInRam(),
					super.getTmpDirectories()
					);
			this.variants.setDestructiveIteration(true);
			
			
			final VCFHeader header=(VCFHeader)r.getHeader();
			
	       SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
			while(r.hasNext())
				{
				final VariantContext ctx= progress.watch(r.next());
				/* discard non SNV variant */
				if(!ctx.isSNP() ||  !ctx.isVariant())
					{
					continue;
					}
				
				final Collection<KnownGene> genes= this.knownGenes.getOverlapping(
						new Interval(ctx.getContig(),ctx.getStart(),ctx.getEnd())
						); 
				for(final KnownGene kg:genes) {
					for(final Allele alt: ctx.getAlternateAlleles()) {
						challenge(ctx,alt,kg);
						}
					}
				}
			progress.finish();
			this.variants.doneAdding();
			final GeneticCode geneticCode= GeneticCode.getStandard();
			
			varIter = this.variants.iterator();
			ArrayList<Variant> buffer= new ArrayList<>();
			for(;;)
				{
				Variant variant = null;
				if(varIter.hasNext())
					{
					variant = varIter.next();
					}
				if(variant==null || !(!buffer.isEmpty() && buffer.get(0).contig.equals(variant.contig) &&  buffer.get(0).transcriptName.equals(variant.transcriptName)))
					{
					if(!buffer.isEmpty()) {
					for(int i=0;i< buffer.size();++i)
						{
						final Variant v1  = buffer.get(i);
						for(int j=i+1;j< buffer.size();++j)
							{
							final Variant v2  = buffer.get(j);
							if(v1.codonStart() != v2.codonStart()) continue;
							if(v1.positionInCodon() == v2.positionInCodon()) continue;
							if(!v1.wildCodon.equals(v2.wildCodon))
								{
								throw new IllegalStateException();
								}
							
							StringBuilder combinedCodon = new StringBuilder(v1.wildCodon);
							combinedCodon.setCharAt(v1.positionInCodon(), v1.mutCodon.charAt(v1.positionInCodon()));
							combinedCodon.setCharAt(v2.positionInCodon(), v2.mutCodon.charAt(v2.positionInCodon()));
							
							final String pwild = new ProteinCharSequence(geneticCode, v1.wildCodon).getString();
							final String p1 = new ProteinCharSequence(geneticCode, v1.mutCodon).getString();
							final String p2 = new ProteinCharSequence(geneticCode, v2.mutCodon).getString();
							final String pCombined = new ProteinCharSequence(geneticCode, combinedCodon).getString();
							if(!pCombined.equals(pwild) &&
								p1.equals(pwild) &&
								p2.equals(pwild)
								) {
								System.out.print(v1.toString());
								System.out.print("\t");
								System.out.print(p1);
								System.out.print("\t");
								System.out.print(v2.toString());
								System.out.print("\t");
								System.out.print(p2);
								System.out.print("\t");
								System.out.print(pCombined);
								System.out.println();


							}
							
							}
						}
					}
					buffer.clear();
					if(variant==null) break;
					}
				buffer.add(variant);
				}
			
			
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			CloserUtil.close(varIter);
			if(this.variants!=null) this.variants.cleanup();
			this.variants=null;
			}
		
		}
	
	private void challenge(
			final VariantContext ctx,
			final Allele allele,
			final KnownGene gene
			) throws IOException
		{
		if(genomicSequence==null || !genomicSequence.getChrom().equals(ctx.getContig()))
			{
			LOG.info("getting genomic Sequence for "+gene.getContig());
			genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, gene.getContig());
			}
				
		Variant variant=null;
    		
		final StringBuilder wildRNA=new StringBuilder();
		final MutedSequence mutRNA=new MutedSequence(wildRNA);
    	
		if(gene.isPositiveStrand())
			{
	    	for(int exon_index=0;exon_index< gene.getExonCount();++exon_index)
	    		{
	    		final KnownGene.Exon exon= gene.getExon(exon_index);
	    		for(int i= exon.getStart();
						i< exon.getEnd();
						++i)
					{
					if(i< gene.getCdsStart()) continue;
					if(i>=gene.getCdsEnd()) break;

					wildRNA.append(genomicSequence.charAt(i));
					
					if(variant==null && ctx.getStart()-1==i)
						{
						variant = new Variant(ctx,allele,gene);
						variant.sorting_id = ID_GENERATOR++;
						variant.position_in_cdna=wildRNA.length()-1;
						char mutBase= allele.getBaseString().charAt(0);
						mutRNA.put( wildRNA.length()-1, mutBase 	);
						}
					}
	    		}
			}
		else
			{
			int exon_index = gene.getExonCount()-1;
			while(exon_index >=0)
				{
				final KnownGene.Exon exon= gene.getExon(exon_index);
				
				for(int i= exon.getEnd()-1;
					    i>= exon.getStart();
					--i)
					{
					if(i>= gene.getCdsEnd()) continue;
					if(i<  gene.getCdsStart()) break;
	
					wildRNA.append(AcidNucleics.complement(genomicSequence.charAt(i)));
					
					if(variant==null && ctx.getStart()-1==i)
						{
						variant = new Variant(ctx,allele,gene);
						variant.sorting_id = ID_GENERATOR++;
						variant.position_in_cdna=wildRNA.length()-1;
						char mutBase=AcidNucleics.complement(allele.getBaseString().charAt(0));
						mutRNA.put( wildRNA.length()-1, mutBase 	);
						}
    				}
    			--exon_index;
    			}
			}
        	
    	if(variant!=null)
    		{
    		variant.wildCodon="";
    		variant.mutCodon="";
    		for(int i=0;i< 3;++i)
    			{
    			int pos = variant.codonStart()+i;
    			variant.wildCodon += (pos< wildRNA.length()?wildRNA.charAt(pos):'*');
    			variant.mutCodon +=  (pos< mutRNA.length()?mutRNA.charAt(pos):'*');
    			}
    		variant.wildCodon = variant.wildCodon.toUpperCase();
    		variant.mutCodon = variant.mutCodon.toUpperCase();
    		
    		if(variant.wildCodon.equals(variant.mutCodon)) {
    			LOG.info("Uh???????");
    			return;
    			}
    		this.variants.add(variant);
    		}
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		if(this.referenceFile==null)
			{
			return wrapException("Undefined REFERENCE. option: -"+OPTION_REFERENCEFILE);
			}
		if(this.kgURI==null || this.kgURI.trim().isEmpty())
			{
			return wrapException("Undefined kgURI. option: -"+OPTION_KGURI);
			}
		VcfIterator iter=null;
		try {
			iter = super.openVcfIterator(inputName);
			return scan(iter);
		} catch (Exception e) {
			return wrapException(e);
			} finally
			{
			CloserUtil.close(iter);	
			}
		}
	
	public static void main(String[] args)
		{
		new VCFStopCodon().instanceMainWithExit(args);
		}
	}