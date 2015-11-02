/**
 * Author:
 * 	Pierre Lindenbaum PhD
 * Date:
 * 	Dec-2010
 * Contact:
 * 	plindenbaum@yahoo.fr
 * Reference:
 *   http://plindenbaum.blogspot.com/2011/01/my-tool-to-annotate-vcf-files.html
 * Motivation:
 * 	Annotate a VCF file with a knownGene structure
 * Compilation:
 */
package com.github.lindenb.jvarkit.tools.vcfannot;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Pattern;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamSequenceRecordTreeMap;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;





/**
 * VCFAnnotator
 * Annotator for VCF
 *
 */
public class VCFAnnotator extends AbstractVCFAnnotator
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(VCFAnnotator.class);
	public static final String TAG="PRED";
	public static enum FORMAT1{TRANSCRIPT,CDNAPOS,PROTPOS,CODON,AA,SEQONTOLOGY};

	@Override
	public Command createCommand()
		{
		return new MyCommand();
		}
	
	static private class MyCommand extends AbstractVCFAnnotator.AbstractVCFAnnotatorCommand
		{


	private SamSequenceRecordTreeMap<KnownGene> knownGenes=null;
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	
	
	
	private static class MutedSequence extends DelegateCharSequence
		{
		private Map<Integer, Character> pos2char=new TreeMap<Integer, Character>();
		MutedSequence(CharSequence wild)
			{
			super(wild);
			}
		
		void put(int pos,char c)
			{
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
		
	
	class Annotation
		{
		KnownGene kg;
		Allele alt2;
		String exon_name="";
		String intron_name="";
		Integer position_cdna=null;
		String wildAA="";
		String mutAA="";
		String wildCodon="";
		String mutCodon="";
		Integer position_protein=null;
		Set<SequenceOntologyTree.Term> seqont=new HashSet<SequenceOntologyTree.Term>();
		
		
		public String toString()
			{
			StringBuilder b=new StringBuilder();
			
			b.append(kg==null?"":kg.getName());
			b.append('|');
			b.append(position_cdna==null?"":String.valueOf(position_cdna));
			b.append('|');
			b.append(position_protein==null?"":String.valueOf(position_protein));
			b.append('|');
			b.append(wildCodon.isEmpty() && mutCodon.isEmpty()?"":wildCodon+"/"+mutCodon);
			b.append('|');
			b.append(wildAA.isEmpty() && mutAA.isEmpty()?"":wildAA+"/"+mutAA);
			b.append("|");
			boolean first=true;
			for(SequenceOntologyTree.Term t:seqont)
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
			
			return b.toString();
			}
		
		}	
	
	private void loadKnownGenesFromUri() throws IOException
		{
		if(this.indexedFastaSequenceFile.getSequenceDictionary()==null)
			{
			throw new IOException("Cannot get sequence dictionary for REF : "+getMessageBundle("picard.dictionary.needed"));
			}
		int n_genes=0;
		this.knownGenes=new SamSequenceRecordTreeMap<KnownGene>(
				this.indexedFastaSequenceFile.getSequenceDictionary()
				);
		LOG.info("loading genes");
		Set<String> unknown=new HashSet<String>();
		BufferedReader in=IOUtils.openURIForBufferedReading(this.kgURI);
		String line;
		Pattern tab=Pattern.compile("[\t]");
		while((line=in.readLine())!=null)
			{
			if(line.isEmpty()) continue;
			String tokens[]=tab.split(line);
			KnownGene g=new KnownGene(tokens);

			if(!this.knownGenes.put(g.getContig(), g.getTxStart()+1, g.getTxEnd(),g))
				{
				if(!unknown.contains(g.getContig()))
					{
					LOG.warn("The reference "+REF+" doesn't contain chromosome "+g.getContig());
					unknown.add(g.getContig());
					}
				}
			else
				{
				++n_genes;
				}
			}
		in.close();
		LOG.info("genes:"+n_genes);
		}
	private boolean isStop(char c)
		{
		return !Character.isLetter(c);
		}
	
	private  boolean isSimpleBase(Allele a)
		{
		String s=a.getBaseString();
		if(s.length()!=1) return false;
		switch(s.charAt(0))
			{
			case 'A':case 'a':case 'T':case 't':
			case 'G':case 'g':case 'C':case 'c':
				return true;
			}
		return false;
		}
	
	
	@Override
		protected Collection<Throwable> doVcfToVcf(String inputName,
				VcfIterator r, VariantContextWriter w) throws IOException {
		GenomicSequence genomicSequence=null;
		LOG.info("opening REF:"+REF);
		this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(REF);
		loadKnownGenesFromUri();
		VCFHeader header=(VCFHeader)r.getHeader();
		
		
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		addMetaData(h2);
		
		StringBuilder format=new StringBuilder();
		for(FORMAT1 f:FORMAT1.values())
			{
			if(format.length()>0) format.append("|"); 
			 format.append(f.name()); 
			}
		
		h2.addMetaDataLine(new VCFInfoHeaderLine(TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,
				"Prediction from "+getClass().getSimpleName()+
				". Format: "+format
				));
		
		
        w.writeHeader(h2);

		
		SequenceOntologyTree soTree=SequenceOntologyTree.getInstance();
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

		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
		while(r.hasNext())
			{
			VariantContext ctx=r.next();
			
			progress.watch(ctx.getContig(), ctx.getStart());
			
			List<KnownGene> genes=this.knownGenes.getOverlapping(
					ctx.getContig(), ctx.getStart(), ctx.getEnd() //1-based
					);
			List<Annotation> ctx_annotations=new ArrayList<Annotation>();
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
				
				for(KnownGene gene:genes)
					{
					GeneticCode geneticCode=GeneticCode.getStandard();
            		
            		
					for(Allele alt2:ctx.getAlternateAlleles())
						{
						Annotation annotations=new Annotation();
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
		        		int position_in_cdna=-1;
		        		
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
		            		if(position < gene.getCdsStart())
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
			            			KnownGene.Exon exon= gene.getExon(exon_index);
			            			
			            			
			            			
			            			
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
			        						position_in_cdna=wildRNA.length();
			        						annotations.position_cdna=position_in_cdna;
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
			            							position_in_cdna,
			            							alt2.getBaseString().charAt(0)
			            							);
			            					
			            					}
			            				
			            				if(wildRNA.length()%3==0 && wildRNA.length()>0 && wildProt==null)
				            				{
				            				wildProt=new ProteinCharSequence(geneticCode,wildRNA);
				            				mutProt=new ProteinCharSequence(geneticCode,mutRNA);
				            				}
			            				}
			            			KnownGene.Intron intron= exon.getNextIntron();
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

		            		if(position < gene.getCdsStart())
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

			            			KnownGene.Exon exon= gene.getExon(exon_index);
			            			
			            			
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
			            					position_in_cdna=wildRNA.length();
			        						annotations.position_cdna=position_in_cdna;
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
				        								position_in_cdna,
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
			            			KnownGene.Intron intron= exon.getPrevIntron();
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
		        			position_in_cdna>=0)
			    			{
		            		int pos_aa=position_in_cdna/3;
		            		int mod= position_in_cdna%3;
		            		annotations.wildCodon=(""+
		            			wildRNA.charAt(position_in_cdna-mod+0)+
		            			wildRNA.charAt(position_in_cdna-mod+1)+
		            			wildRNA.charAt(position_in_cdna-mod+2)
		            			);
		            		annotations.mutCodon=(""+
		            			mutRNA.charAt(position_in_cdna-mod+0)+
		            			mutRNA.charAt(position_in_cdna-mod+1)+
		            			mutRNA.charAt(position_in_cdna-mod+2)
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
			
		
			
			Set<String> info=new HashSet<String>(ctx_annotations.size());
			for(Annotation a:ctx_annotations)
				{
				info.add(a.toString());
				}
			
			VariantContextBuilder vb=new VariantContextBuilder(ctx);
			vb.attribute(TAG, info.toArray());
			w.add(vb.make());
			}
		CloserUtil.close(this.indexedFastaSequenceFile);
		this.indexedFastaSequenceFile=null;
		return RETURN_OK;
		}
	
		
		@Override
			protected Collection<Throwable> call(String inputName) throws Exception
			{
			if(this.REF==null)
				{
				return wrapException("Undefined REFERENCE.");
				}
			
			return doVcfToVcf(inputName);
			}
		}
	
	public static void main(String[] args)
		{
		new VCFAnnotator().instanceMainWithExit(args);
		}
	}