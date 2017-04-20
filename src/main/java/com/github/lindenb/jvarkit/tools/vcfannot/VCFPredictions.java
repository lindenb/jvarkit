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
package com.github.lindenb.jvarkit.tools.vcfannot;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Pattern;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


/**
 * VCFAnnotator
 * Annotator for VCF
 *
 */
@Program(name="vcfpredictions",description="Basic Variant Effect prediction using ucsc-known gene")
public class VCFPredictions extends Launcher
	{
	private static final Logger LOG = Logger.build(VCFPredictions.class).make();

	private IntervalTreeMap<List<KnownGene>> knownGenes=null;
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	



	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;


	@Parameter(names={"-k","--knownGene"},description="KnownGene data URI/File. Beware chromosome names are formatted the same as your REFERENCE.",required=true)
	private String kgURI = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz";

	@Parameter(names={"-soacn","--printsoacn"},description="Print SO:term accession rather than label")
	private boolean print_SO_ACN = false;

	@Parameter(names={"-vep","--vep"},description="Variant Effect Predictor output Syntax")
	private boolean vepSyntax = false;

	@Parameter(names={"-R","--reference"},description="Indexed fasta Reference",required=true)
	private File referenceFile = null;

	
	private static class MutedSequence extends DelegateCharSequence
		{
		private final Map<Integer, Character> pos2char=new TreeMap<Integer, Character>();
		MutedSequence(final CharSequence wild)
			{
			super(wild);
			}
		
		void put(final int pos,final char c)
			{
			this.pos2char.put(pos, c);
			}
		
		@Override
		public char charAt(final int i)
			{
			Character c= pos2char.get(i);
			return c==null?getDelegate().charAt(i):c;
			}
		
		}
	
	
	private static class ProteinCharSequence extends DelegateCharSequence
		{
		private GeneticCode geneticCode;
		ProteinCharSequence(final GeneticCode geneticCode,final CharSequence cDNA)
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
		
	
	class Annotation
		{
		KnownGene kg;
		Allele alt2;
		String exon_name="";
		String intron_name="";
		Integer position_cds=null;
		String wildAA="";
		String mutAA="";
		String wildCodon="";
		String mutCodon="";
		Integer position_protein=null;
		Set<SequenceOntologyTree.Term> seqont=new HashSet<SequenceOntologyTree.Term>();
		
		
		public String toString()
			{
			final StringBuilder b=new StringBuilder();
			if(VCFPredictions.this.vepSyntax) {
				boolean first=true;
				//Allele|Feature|Feature_type|Consequence|CDS_position|Protein_position|Amino_acids|Codons
				b.append(alt2==null?"":alt2.getBaseString());
				b.append('|');
				b.append(kg==null?"":kg.getName());
				b.append('|');
				b.append("Transcript");
				b.append('|');
				for(final SequenceOntologyTree.Term t:seqont)
					{
					if(!first) b.append('&');
					first=false;
					b.append(t.getLabel().replaceAll("[ ]","_"));
						
					}
				b.append('|');
				b.append(position_cds==null?"":String.valueOf(position_cds+1));
				b.append('|');
				b.append(position_protein==null?"":String.valueOf(position_protein));
				b.append('|');
				b.append(wildAA.isEmpty() && mutAA.isEmpty()?"":wildAA+"/"+mutAA);
				b.append('|');
				b.append(wildCodon.isEmpty() && mutCodon.isEmpty()?"":wildCodon+"/"+mutCodon);
				}
			else {
			b.append(kg==null?"":kg.getName());
			b.append('|');
			b.append(position_cds==null?"":String.valueOf(position_cds));
			b.append('|');
			b.append(position_protein==null?"":String.valueOf(position_protein));
			b.append('|');
			b.append(wildCodon.isEmpty() && mutCodon.isEmpty()?"":wildCodon+"/"+mutCodon);
			b.append('|');
			b.append(wildAA.isEmpty() && mutAA.isEmpty()?"":wildAA+"/"+mutAA);
			b.append("|");
			boolean first=true;
			for(final SequenceOntologyTree.Term t:seqont)
				{
				if(!first) b.append('&');
				first=false;
				if(print_SO_ACN)
					{
					b.append(t.getAcn());
					}
				else
					{
					b.append(t.getLabel().replaceAll("[ ]","_"));
					}
				}
			}
			return b.toString();
			}
		
		}	
	
	private void loadKnownGenesFromUri() throws IOException
		{
		BufferedReader in=null;
		try {
			if (this.indexedFastaSequenceFile.getSequenceDictionary() == null) {
				throw new IOException(
						"Cannot get sequence dictionary for REF : " + getMessageBundle("picard.dictionary.needed"));
			}
			int n_genes = 0;
			this.knownGenes = new IntervalTreeMap<>();
			LOG.info("loading genes");
			in = IOUtils.openURIForBufferedReading(this.kgURI);
			String line;
			final Pattern tab = Pattern.compile("[\t]");
			while ((line = in.readLine()) != null) {
				if (line.isEmpty())
					continue;
				final String tokens[] = tab.split(line);
				final KnownGene g = new KnownGene(tokens);
				if (this.indexedFastaSequenceFile.getSequenceDictionary().getSequence(g.getContig()) == null) {
					continue;
				}
				final int extend_gene_search = 5000; // because we want to set
														// SO:5KB_upstream_variant

				final Interval interval = new Interval(g.getContig(),
						Math.max(1, g.getTxStart() + 1 - extend_gene_search), g.getTxEnd() + extend_gene_search);
				List<KnownGene> L= this.knownGenes.get(interval);
				if(L==null) {
					L=new ArrayList<>(2);
					this.knownGenes.put(interval, L);
				}
				
				L.add(g);
			}
			in.close();
			in = null;
			LOG.info("genes:" + n_genes);
			}
		finally {
			CloserUtil.close(in);
			}
		}
	private boolean isStop(char c)
		{
		return !Character.isLetter(c);
		}
	
	private  boolean isSimpleBase(final Allele a)
		{
		if(a.isSymbolic()) return false;
		final String s=a.getBaseString();
		if(s.length()!=1) return false;
		switch(s.charAt(0))
			{
			case 'A':case 'a':case 'T':case 't':
			case 'G':case 'g':case 'C':case 'c':
				return true;
			}
		return false;
		}
	
	public static final String TAG="PRED";
	public static enum FORMAT1{TRANSCRIPT,CDSPOS,PROTPOS,CODON,AA,SEQONTOLOGY};
	
	
	@Override
	protected int doVcfToVcf(final String inputName, final VcfIterator r, VariantContextWriter w)
			 {
		GenomicSequence genomicSequence=null;
		try {
		LOG.info("opening REF:"+this.referenceFile);
		this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(this.referenceFile);
		loadKnownGenesFromUri();
		final VCFHeader header=(VCFHeader)r.getHeader();
		
		
		final VCFHeader h2=new VCFHeader(header);
		addMetaData(h2);
		
		
		if(this.vepSyntax)
			{
			h2.addMetaDataLine(new VCFInfoHeaderLine("CSQ",
					VCFHeaderLineCount.UNBOUNDED,
					VCFHeaderLineType.String,
					"Consequence type as predicted by VEP"+
					". Format: Allele|Feature|Feature_type|Consequence|CDS_position|Protein_position|Amino_acids|Codons"
					));
			}
		else
			{
			final StringBuilder format=new StringBuilder();
			for(FORMAT1 f:FORMAT1.values())
				{
				if(format.length()>0) format.append("|"); 
				 format.append(f.name()); 
				}
			
			h2.addMetaDataLine(new VCFInfoHeaderLine(TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,
					"Prediction from "+getClass().getSimpleName()+
					". Format: "+format
					));
			}
		
        w.writeHeader(h2);

		
		final SequenceOntologyTree soTree=SequenceOntologyTree.getInstance();
		final SequenceOntologyTree.Term so_intron=soTree.getTermByAcn("SO:0001627");
		final SequenceOntologyTree.Term so_exon=soTree.getTermByAcn("SO:0001791");
		final SequenceOntologyTree.Term so_splice_donor=soTree.getTermByAcn("SO:0001575");
		final SequenceOntologyTree.Term so_splice_acceptor=soTree.getTermByAcn("SO:0001574");
		final SequenceOntologyTree.Term so_5_prime_UTR_variant=soTree.getTermByAcn("SO:0001623");
		final SequenceOntologyTree.Term so_3_prime_UTR_variant=soTree.getTermByAcn("SO:0001624");
		final SequenceOntologyTree.Term so_splicing_variant=soTree.getTermByAcn("SO:0001568");
		final SequenceOntologyTree.Term so_stop_lost=soTree.getTermByAcn("SO:0001578");
		final SequenceOntologyTree.Term so_stop_gained=soTree.getTermByAcn("SO:0001587");
		final SequenceOntologyTree.Term so_coding_synonymous=soTree.getTermByAcn("SO:0001819");
		final SequenceOntologyTree.Term so_coding_non_synonymous=soTree.getTermByAcn("SO:0001583");
		final SequenceOntologyTree.Term so_intergenic=soTree.getTermByAcn("SO:0001628");
		final SequenceOntologyTree.Term so_nc_transcript_variant=soTree.getTermByAcn("SO:0001619");
		final SequenceOntologyTree.Term so_non_coding_exon_variant=soTree.getTermByAcn("SO:0001792");
		final SequenceOntologyTree.Term _2KB_upstream_variant=soTree.getTermByAcn("SO:0001636");
		final SequenceOntologyTree.Term _5KB_upstream_variant=soTree.getTermByAcn("SO:0001635");
		final SequenceOntologyTree.Term _5KB_downstream_variant=soTree.getTermByAcn("SO:0001633");
		final SequenceOntologyTree.Term _500bp_downstream_variant=soTree.getTermByAcn("SO:0001634");
		
		
		final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
		while(r.hasNext())
			{
			final VariantContext ctx=progress.watch(r.next());
			
			final List<KnownGene> genes=new ArrayList<>();
			
			for(final List<KnownGene> l2: this.knownGenes.getOverlapping(new Interval(
					ctx.getContig(),
					ctx.getStart(),
					ctx.getEnd() //1-based
					)))
				{
				genes.addAll(l2);
				}
			final List<Annotation> ctx_annotations=new ArrayList<Annotation>();
			if(genes==null || genes.isEmpty())
				{
				//intergenic
				Annotation a=new Annotation();
				a.seqont.add(so_intergenic);
				ctx_annotations.add(a);
				}
			else
				{
				if(genomicSequence==null || !genomicSequence.getChrom().equals(ctx.getContig()))
					{
					LOG.info("getting genomic Sequence for "+ctx.getContig());
					genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, ctx.getContig());
					}
				
				for(final KnownGene gene:genes)
					{
					final GeneticCode geneticCode=GeneticCode.getStandard();
            		
            		
					for(final Allele alt2:ctx.getAlternateAlleles())
						{
						if(alt2.isSymbolic()) continue;
						if(alt2.isReference()) continue;
						final Annotation annotations=new Annotation();
						annotations.kg=gene;
						annotations.alt2=alt2;
						
						if(gene.isNonCoding())
							{
							annotations.seqont.add(so_nc_transcript_variant);
							continue;
							}
						
						
						ctx_annotations.add(annotations);

		        		StringBuilder wildRNA=null;
		        		ProteinCharSequence wildProt=null;
		        		ProteinCharSequence mutProt=null;
		        		MutedSequence mutRNA=null;
		        		int position_in_cds=-1;
		        		
		        		final int position=ctx.getStart()-1;
		        		if(!String.valueOf(genomicSequence.charAt(position)).equalsIgnoreCase(ctx.getReference().getBaseString()))
		        			{
		        			if(isSimpleBase(ctx.getReference()))
			        			{
			        			LOG.warn("Warning REF!=GENOMIC SEQ!!! at "+position+"/"+ctx.getReference());
			        			}
		        			continue;
		        			}
		        		
		        		if(gene.isPositiveStrand())
		            		{
		        			if(position < gene.getTxStart() - 2000) {
		        				annotations.seqont.add(_5KB_upstream_variant);
		        				}
		        			else if(position < gene.getTxStart()) {
		        				annotations.seqont.add(_2KB_upstream_variant);
		        				}
		        			else if( position >= gene.getTxEnd() + 500) {
		        				annotations.seqont.add(_5KB_downstream_variant);
		        				}
		        			else if( position >= gene.getTxEnd() ) {
		        				annotations.seqont.add(_500bp_downstream_variant);
		        				}
		        			else if(position < gene.getCdsStart())
		            			{
		            			annotations.seqont.add(so_5_prime_UTR_variant);//UTR5
		            			}
		            		else if( gene.getCdsEnd()<= position )
		            			{
		            			annotations.seqont.add(so_3_prime_UTR_variant);
		            			}
		            		else
			            		{
			            		int exon_index=0;
			            		while(exon_index< gene.getExonCount())
			            			{
			            			final KnownGene.Exon exon= gene.getExon(exon_index);
			            			
			            			for(int i= exon.getStart();
			            					i< exon.getEnd();
			            					++i)
			            				{
			            				
			            				if(i==position)
			        						{
			        						annotations.exon_name= exon.getName();
			        						if(exon.isNonCoding())
			            						{
			            						annotations.seqont.add(so_non_coding_exon_variant);
			            						}
			        						}
			            				if(i< gene.getTxStart()) continue;
			            				if(i< gene.getCdsStart()) continue;
			            				if(i>=gene.getCdsEnd()) break;
			        					
			        					if(wildRNA==null)
			        						{
			        						wildRNA=new StringBuilder();
			        						mutRNA=new MutedSequence(wildRNA);
			        						}
			        					
			        					if(i==position)
			        						{
			        						annotations.seqont.add(so_exon);
			        						annotations.exon_name=exon.getName();
			        						position_in_cds=wildRNA.length();
			        						annotations.position_cds= position_in_cds;
			        						//in splicing ?
			        						if(exon.isSplicing(position))
			        							{
			        							
			        							if(exon.isSplicingAcceptor(position))
			        								{
			        								annotations.seqont.add(so_splice_acceptor); //SPLICING_ACCEPTOR
			        								}
			        							else  if(exon.isSplicingDonor(position))
			        								{
			        								annotations.seqont.add(so_splice_donor); // SPLICING_DONOR
			        								}
			        							else //??
			        								{
				        							annotations.seqont.add(so_splicing_variant);
			        								}
			        							}
			        						}
			        					
			            				wildRNA.append(genomicSequence.charAt(i));
			            				
			            				if(i==position && 
			            						isSimpleBase(alt2) && 
			            						isSimpleBase(ctx.getReference()))
			            					{
			            					mutRNA.put(
			            							position_in_cds,
			            							alt2.getBaseString().charAt(0)
			            							);
			            					
			            					}
			            				
			            				if(wildRNA.length()%3==0 && wildRNA.length()>0 && wildProt==null)
				            				{
				            				wildProt=new ProteinCharSequence(geneticCode,wildRNA);
				            				mutProt=new ProteinCharSequence(geneticCode,mutRNA);
				            				}
			            				}
			            			final KnownGene.Intron intron= exon.getNextIntron();
			            			if(intron!=null && intron.contains(position))
			            				{
			            				annotations.intron_name=intron.getName();
			            				annotations.seqont.add(so_intron);
			            				
			            				if(intron.isSplicing(position))
			        						{
			        						if(intron.isSplicingAcceptor(position))
			        							{
			        							annotations.seqont.add(so_splice_acceptor);
			        							}
			        						else if(intron.isSplicingDonor(position))
			        							{
			        							annotations.seqont.add(so_splice_donor);
			        							}
			        						else //???
			        							{
			        							annotations.seqont.add(so_splicing_variant);
			        							}
			        						}
			            				}
			            			++exon_index;
			            			}
			            		}
		            		
		            		
		            		}
		            	else // reverse orientation
		            		{
		            		if(position >= gene.getTxEnd() + 2000) {
		        				annotations.seqont.add(_5KB_upstream_variant);
		        				}
		        			else if(position >= gene.getTxEnd()) {
		        				annotations.seqont.add(_2KB_upstream_variant);
		        				}
		        			else if( position < gene.getTxStart() - 500) {
		        				annotations.seqont.add(_5KB_downstream_variant);
		        				}
		        			else if( position < gene.getTxStart() ) {
		        				annotations.seqont.add(_500bp_downstream_variant);
		        				}
		        			else if(position < gene.getCdsStart())
		            			{
		            			annotations.seqont.add(so_3_prime_UTR_variant);
		            			}
		            		else if( gene.getCdsEnd()<=position )
		            			{
		            			annotations.seqont.add(so_5_prime_UTR_variant);
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
			            				
			            				
			            				if(i==position)
			        						{
			            					annotations.exon_name=exon.getName();
			            					if(exon.isNonCoding())
			            						{
			            						annotations.seqont.add(so_non_coding_exon_variant);
			            						}
			        						}
			            				if(i>= gene.getCdsEnd()) continue;
			            				if(i<  gene.getCdsStart()) break;
			            				
			            				
			            				if(wildRNA==null)
			        						{
			        						wildRNA=new StringBuilder();
			        						mutRNA=new MutedSequence(wildRNA);
			        						}
			            				
			            				if(i==position)
			        						{
			            					annotations.seqont.add(so_exon);
			            					position_in_cds=wildRNA.length();
			        						annotations.position_cds=position_in_cds;
			        						//in splicing ?
			        						if(exon.isSplicing(position))
			        							{		        							
			        							if(exon.isSplicingAcceptor(position))
			        								{
			        								annotations.seqont.add(so_splice_acceptor);
			        								}
			        							else  if(exon.isSplicingDonor(position))
			        								{
			        								annotations.seqont.add(so_splice_donor);
			        								}
			        							else //?
			        								{
			        								annotations.seqont.add(so_splicing_variant);
			        								}
			        							}
			        						
			        						if(isSimpleBase(alt2) &&
			        							isSimpleBase(ctx.getReference()))
				        						{
				        						mutRNA.put(
				        								position_in_cds,
				        								AcidNucleics.complement(alt2.getBaseString().charAt(0))
				        								);
				        						}
			        						}
			            				
			            				wildRNA.append(AcidNucleics.complement(genomicSequence.charAt(i)));
			            				if( wildRNA.length()%3==0 &&
			            					wildRNA.length()>0 &&
			            					wildProt==null)
				            				{
				            				wildProt=new ProteinCharSequence(geneticCode,wildRNA);
				            				mutProt=new ProteinCharSequence(geneticCode,mutRNA);
				            				}
			            				}
			            			final KnownGene.Intron intron= exon.getPrevIntron();
			            			if(intron!=null &&
			            				intron.contains(position))
			            				{
			            				annotations.intron_name=intron.getName();
			            				annotations.seqont.add(so_intron);
			            				
			            				if(intron.isSplicing(position))
			        						{
			        						if(intron.isSplicingAcceptor(position))
			        							{
			        							annotations.seqont.add(so_splice_acceptor);
			        							}
			        						else if(intron.isSplicingDonor(position))
			        							{
			        							annotations.seqont.add(so_splice_donor);
			        							}
			        						else //?	
			        							{
			        							annotations.seqont.add(so_splicing_variant);
			        							}
			        						}
			            				}
			            			--exon_index;
			            			}
			            		}

		            		}//end of if reverse
		        		
		        		
		        		if( isSimpleBase(alt2) &&
		        			isSimpleBase(ctx.getReference()) &&
		        			wildProt!=null &&
		        			mutProt!=null && 
		        			position_in_cds>=0)
			    			{
		            		final int pos_aa=position_in_cds/3;
		            		final int mod= position_in_cds%3;
		            		annotations.wildCodon=(""+
		            			wildRNA.charAt(position_in_cds-mod+0)+
		            			wildRNA.charAt(position_in_cds-mod+1)+
		            			wildRNA.charAt(position_in_cds-mod+2)
		            			);
		            		annotations.mutCodon=(""+
		            			mutRNA.charAt(position_in_cds-mod+0)+
		            			mutRNA.charAt(position_in_cds-mod+1)+
		            			mutRNA.charAt(position_in_cds-mod+2)
		            			);
		            		annotations.position_protein=(pos_aa+1);
		            		annotations.wildAA=String.valueOf(wildProt.charAt(pos_aa));
		            		annotations.mutAA=(String.valueOf(mutProt.charAt(pos_aa)));
		            		
		            		annotations.seqont.remove(so_exon);
		            		
			    			if(isStop(wildProt.charAt(pos_aa)) &&
			    			   !isStop(mutProt.charAt(pos_aa)))
			    				{
			    				annotations.seqont.add(so_stop_lost);
			    				}
			    			else if( !isStop(wildProt.charAt(pos_aa)) &&
			    				 isStop(mutProt.charAt(pos_aa)))
			    				{
			    				annotations.seqont.add(so_stop_gained);
			    				}
			    			else if(wildProt.charAt(pos_aa)==mutProt.charAt(pos_aa))
			    				{
			    				annotations.seqont.add(so_coding_synonymous);
			    				}
			    			else
			    				{
			    				annotations.seqont.add(so_coding_non_synonymous);
			    				}
			    			}
		        		
						}
					}
				}
			
		
			
			final Set<String> info=new HashSet<String>(ctx_annotations.size());
			for(final Annotation a:ctx_annotations)
				{
				info.add(a.toString());
				}
			
			final VariantContextBuilder vb=new VariantContextBuilder(ctx);
			vb.attribute((this.vepSyntax?"CSQ":TAG), info.toArray());
			w.add(vb.make());
			}
		
		return RETURN_OK;
		} catch(Exception err ) {
			return wrapException(err);
		} finally {
			CloserUtil.close(this.indexedFastaSequenceFile);
		}
		}
	
	@Override
	public int doWork(final List<String> args) {
			if(this.referenceFile==null) 
			{
			return wrapException("Reference undefined.");
			}
		if(this.kgURI==null || this.kgURI.trim().isEmpty()) 
			{
			return wrapException("knownGene undefined.");
			}
		return doVcfToVcf(args,outputFile);
		}
	
	public static void main(String[] args)
		{
		new VCFPredictions().instanceMainWithExit(args);
		}
	}