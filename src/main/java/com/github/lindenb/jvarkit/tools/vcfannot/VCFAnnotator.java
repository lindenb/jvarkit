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
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Pattern;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMSequenceRecord;

import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.IOUtils;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter;



class MutedSequence extends DelegateCharSequence
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


class ProteinCharSequence extends DelegateCharSequence
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


/**
 * VCFAnnotator
 * Annotator for VCF
 *
 */
public class VCFAnnotator extends AbstractVCFFilter
	{
	static Log LOG=Log.getInstance(VCFAnnotator.class);
	public File IN;
	public File REF;
	public String KG;
	
	private IntervalTreeMap<KnownGene> knownGenes=new IntervalTreeMap<KnownGene>();
	private IndexedFastaSequenceFile indexedFastaSequenceFile;
	
	
	static class Annotation
		{
		KnownGene kg;
		Allele alt;
		String exon_name="";
		String intron_name="";
		Integer position_cdna=null;
		String wildAA="";
		String mutAA="";
		String wildCodon="";
		String mutCodon="";
		Integer position_protein;
		Set<SequenceOntologyTree.Term> seqont=new HashSet<SequenceOntologyTree.Term>();
		}	
	
	private void loadKnownGenes()throws IOException
		{
		LOG.info("loading genes");
		Set<String> unseen=new HashSet<String>();
		BufferedReader in=IOUtils.openURIForBufferedReading(KG);
		String line;
		Pattern tab=Pattern.compile("[\t]");
		while((line=in.readLine())!=null)
			{
			if(line.isEmpty()) continue;
			String tokens[]=tab.split(line);
			KnownGene g=new KnownGene(tokens);
			SAMSequenceRecord ssr=indexedFastaSequenceFile.getSequenceDictionary().getSequence(g.getChr());
			if(ssr==null)
				{
				if(unseen.add(g.getChr()))
					{
					LOG.error("The reference "+REF+" doesn't contain chromosome "+g.getChr());
					}
				continue;
				}
			Interval interval1=new Interval(g.getChr(), g.getTxStart()+1, g.getTxEnd());
			this.knownGenes.put(interval1, g);
			}
		in.close();
		LOG.info("genes:"+this.knownGenes.size());
		}
	private boolean isStop(char c)
		{
		return !Character.isLetter(c);
		}
	@Override
	protected void doWork(LineReader in, VariantContextWriter w)
		throws IOException
		{	
		LOG.info("opening REF:"+REF);
		this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(REF);
		loadKnownGenes();
		VCFCodec vcfCodec=new VCFCodec();
		VCFHeader header=(VCFHeader)vcfCodec.readHeader(in);
		
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName(),"Filtered with "+getCommandLine()));
		
		
        w.writeHeader(h2);

		
		String line;
		SequenceOntologyTree seqOntologyTree=SequenceOntologyTree.getInstance();
		while((line=in.readLine())!=null)
			{
			VariantContext ctx=vcfCodec.decode(line);
			Collection<KnownGene> genes=this.knownGenes.getOverlapping(
					new Interval(ctx.getChr(), ctx.getStart(), ctx.getEnd()) //1-based
					);
			if(genes==null || genes.isEmpty())
				{
				//continue;
				//TODO if no annotation, create one INTERGENIC
				}
			else
				{
				GenomicSequence genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, ctx.getChr());
				for(KnownGene gene:genes)
					{
            		GeneticCode geneticCode=GeneticCode.getStandard();
            		
            		
					for(Allele alt:ctx.getAlternateAlleles())
						{
						Annotation annotations=new Annotation();
						annotations.kg=gene;
						annotations.alt=alt;
						
						if(!alt.getBaseString().matches("[ATGC]")) continue;
						

		        		StringBuilder wildRNA=null;
		        		ProteinCharSequence wildProt=null;
		        		ProteinCharSequence mutProt=null;
		        		MutedSequence mutRNA=null;
		        		int position_in_cdna=-1;
		        		
		        		final int position=ctx.getStart()-1;
		        		if(!String.valueOf(genomicSequence.charAt(position)).equalsIgnoreCase(ctx.getReference().getBaseString()))
		        			{
		        			LOG.warn("Warning REF!=GENOMIC SEQ!!! at "+genomicSequence+"/"+ctx.getReference());
		        			continue;
		        			}
		        		
		        		if(gene.isPositiveStrand())
		            		{
		            		if(position < gene.getCdsStart())
		            			{
		            			annotations.seqont.add(seqOntologyTree.getTermByAcn("??"));//UTR5 to do
		            			}
		            		else if( gene.getCdsEnd()<= position )
		            			{
		            			annotations.seqont.add(seqOntologyTree.getTermByAcn("??"));//UTR3 to do
		            			}
		            		
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
		        						annotations.seqont.add(null);//TODO put("type", "EXON");
		        						annotations.exon_name=exon.getName();
		        						position_in_cdna=wildRNA.length();
		        						annotations.position_cdna=position_in_cdna;
		        						//in splicing ?
		        						if(exon.isSplicing(position))
		        							{
		        							annotations.seqont.add(null);//TODO put(KEY_SPLICING, "SPLICING");
		        							
		        							if(exon.isSplicingAcceptor(position))
		        								{
		        								annotations.seqont.add(null);//TODO put(KEY_SPLICING, "SPLICING_ACCEPTOR");
		        								}
		        							else  if(exon.isSplicingDonor(position))
		        								{
		        								annotations.seqont.add(null);//TODO put(KEY_SPLICING, "SPLICING_DONOR");
		        								}
		        							}
		        						}
		        					
		            				wildRNA.append(genomicSequence.charAt(i));
		            				
		            				if(i==position)
		            					{
		            					mutRNA.put(
		            							position_in_cdna,
		            							alt.getBaseString().charAt(0)
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
		            				annotations.seqont.add(null);//TODO put("intron.name",intron.getName());
		            				annotations.seqont.add(null);//TODO put("type", "INTRON");
		            				
		            				if(intron.isSplicing(position))
		        						{
		            					annotations.seqont.add(null);//TODO put(KEY_SPLICING, "INTRON_SPLICING");
		        						if(intron.isSplicingAcceptor(position))
		        							{
		        							annotations.seqont.add(null);//TODO put(KEY_SPLICING, "INTRON_SPLICING_ACCEPTOR");
		        							}
		        						else if(intron.isSplicingDonor(position))
		        							{
		        							annotations.seqont.add(null);//TODO put(KEY_SPLICING, "INTRON_SPLICING_DONOR");
		        							}
		        						}
		            				}
		            			++exon_index;
		            			}
		            		
		            		
		            		
		            		}
		            	else // reverse orientation
		            		{
		            	
		            		if(position < gene.getCdsStart())
		            			{
		            			annotations.seqont.add(null);//TODO put(KEY_TYPE, "UTR3");
		            			}
		            		else if( gene.getCdsEnd()<=position )
		            			{
		            			annotations.seqont.add(null);//TODO put(KEY_TYPE, "UTR5");
		            			}
		            	
		            		
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
		            					annotations.seqont.add(null);//TODO put(KEY_TYPE, "EXON");
		            					position_in_cdna=wildRNA.length();
		        						annotations.position_cdna=position_in_cdna;
		        						//in splicing ?
		        						if(exon.isSplicing(position))
		        							{
		        							annotations.seqont.add(null);//TODO put(KEY_SPLICING, "INTRON_SPLICING");
		        							
		        							if(exon.isSplicingAcceptor(position))
		        								{
		        								annotations.seqont.add(null);//TODO put(KEY_SPLICING, "INTRON_SPLICING_ACCEPTOR");
		        								}
		        							else  if(exon.isSplicingDonor(position))
		        								{
		        								annotations.seqont.add(null);//TODO put(KEY_SPLICING, "INTRON_SPLICING_DONOR");
		        								}
		        							}
		        						
		        						
		        						mutRNA.put(
		        								position_in_cdna,
		        								AcidNucleics.complement(alt.getBaseString().charAt(0))
		        								);
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
		            				annotations.seqont.add(null);//TODO .put("intron.name",intron.getName());
		            				annotations.seqont.add(null);//TODO .put(KEY_TYPE, "INTRON");
		            				
		            				if(intron.isSplicing(position))
		        						{
		            					annotations.seqont.add(null);//TODO .put(KEY_SPLICING, "INTRON_SPLICING");
		        						if(intron.isSplicingAcceptor(position))
		        							{
		        							annotations.seqont.add(null);//TODO .put(KEY_SPLICING, "INTRON_SPLICING_ACCEPTOR");
		        							}
		        						else if(intron.isSplicingDonor(position))
		        							{
		        							annotations.seqont.add(null);//TODO .put(KEY_SPLICING, "INTRON_SPLICING_DONOR");
		        							}
		        						}
		            				}
		            			--exon_index;
		            			}

		            		}//end of if reverse
		        		if( wildProt!=null &&
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
			    			if(isStop(wildProt.charAt(pos_aa)) &&
			    			   !isStop(mutProt.charAt(pos_aa)))
			    				{
			    				annotations.seqont.add(null);//TODO .put("type", "EXON_STOP_LOST");
			    				}
			    			else if( !isStop(wildProt.charAt(pos_aa)) &&
			    				 isStop(mutProt.charAt(pos_aa)))
			    				{
			    				annotations.seqont.add(null);//TODO .put("type", "EXON_STOP_GAINED");
			    				}
			    			else if(wildProt.charAt(pos_aa)==mutProt.charAt(pos_aa))
			    				{
			    				annotations.seqont.add(null);//TODO .put("type", "EXON_CODING_SYNONYMOUS");
			    				}
			    			else
			    				{
			    				annotations.seqont.add(null);//TODO .put("type", "EXON_CODING_NON_SYNONYMOUS");
			    				}
			    			}
		        		
						}
					}
				}
			
			VariantContextBuilder vb=new VariantContextBuilder(ctx);
			w.add(vb.make());
			}
			
		}
	
	
	
	public static void main(String[] args)
		{
		new VCFAnnotator().instanceMainWithExit(args);
		}
	}