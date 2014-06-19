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

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Pattern;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;





/**
 * VCFStopCodon
 * @SolenaLS 's idea: consecutive synonymous bases give a stop codon
 *
 */
public class VCFStopCodon extends AbstractVCFFilter2
	{
	private static final String DEFAULT_KG_URI="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz";
	private File REF=null;
	private String kgURI=DEFAULT_KG_URI;

	private Map<String,List<KnownGene>> knownGenes=null;
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private GenomicSequence genomicSequence=null;
	
	
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
		if(this.knownGenes!=null) return;
		int n_genes=0;
		SAMSequenceDictionary dict=this.indexedFastaSequenceFile.getSequenceDictionary();
		this.knownGenes=new HashMap<String, List<KnownGene>>(dict.size());
		for(SAMSequenceRecord ssr:dict.getSequences())
			{
			this.knownGenes.put(ssr.getSequenceName(), new ArrayList<KnownGene>());
			}
		info("loading genes from "+this.kgURI);
		Set<String> unknown=new HashSet<String>();
		LineIterator iter=IOUtils.openURIForLineIterator(this.kgURI);
		Pattern tab=Pattern.compile("[\t]");
		while(iter.hasNext())
			{
			String line=iter.next();
			if(line.isEmpty()) continue;
			String tokens[]=tab.split(line);
			KnownGene g=new KnownGene(tokens);
			if(g.isNonCoding()) continue;

			List<KnownGene> L=this.knownGenes.get(g.getChr());
			if(L==null)
				{
				if(!unknown.contains(g.getChr()))
					{
					warning("The reference "+REF+" doesn't contain chromosome "+g.getChr());
					unknown.add(g.getChr());
					}
				continue;
				}
			L.add(g);
			++n_genes;
			}
		CloserUtil.close(iter);
		for(String C:knownGenes.keySet())
			{
			Collections.sort(knownGenes.get(C),new Comparator<KnownGene>()
				{
				@Override
				public int compare(KnownGene o1, KnownGene o2)
					{
					return o1.getTxStart()-o2.getTxStart();
					}
				});
			}
		info("genes:"+n_genes);
		}
	
	
	
	private static class Variant
		{
		VariantContext ctx;
		int position_in_cdna=-1;
		MutedSequence mutRNA=null;
		Set<String> set=new HashSet<String>();
		@Override
		public String toString() {
			return ctx.getChr()+":"+ctx.getStart()+" "+ctx.getReference()+"/"+ctx.getAlternateAllele(0);
			}
		}
	static final String TAG="STREAMCODON";
	
	protected void dump(VariantContextWriter w,Variant var)
		{
		if(var.set.isEmpty())
			{
			w.add(var.ctx);
			return;
			}
		VariantContextBuilder vcb=new VariantContextBuilder(var.ctx);
		vcb.attribute(TAG, var.set.toArray());
		w.add(vcb.make());
		}
	
	@Override
	protected void doWork(VcfIterator r, VariantContextWriter w)
		throws IOException
		{	
		info("opening REF:"+REF);
	
		this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(REF);
		loadKnownGenesFromUri();
		
		VCFHeader header=(VCFHeader)r.getHeader();
		
		
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",getVersion()));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		h2.addMetaDataLine(new VCFInfoHeaderLine(TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,
				"Prediction from "+getClass().getSimpleName()
				));

		
		
		
        w.writeHeader(h2);
        long n_variants=0L;
        String currChrom=null;
        LinkedList<Variant> variantStack=new LinkedList<Variant>();
        final SAMSequenceDictionary dict=this.indexedFastaSequenceFile.getSequenceDictionary();
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
		List<KnownGene> currListOfGenes=null; 
		while(r.hasNext())
			{
			VariantContext ctx=r.next();
			progress.watch(ctx.getChr(), ctx.getStart());
			
			if(n_variants++%1000==0)
				{
				if(currListOfGenes!=null) info("Genes: "+currListOfGenes.size());
				info("Variants: "+variantStack.size());
				}
			
			if(!ctx.isSNP() || ctx.getAlternateAlleles().size()!=1)
				{
				w.add(ctx);
				continue;
				}

			
			//unknown chromosome
			if(!this.knownGenes.containsKey(ctx.getChr()))
				{
				while(!variantStack.isEmpty())
					{
					dump(w,variantStack.removeFirst());
					}
				warning("unknown chrom "+ctx.getChr());
				currChrom=null;
				continue;
				}
			
			//not same chromosome
			if(!ctx.getChr().equals(currChrom))
				{
				while(!variantStack.isEmpty())
					{
					dump(w,variantStack.removeFirst());
					}
				currChrom=ctx.getChr();
				currListOfGenes=new ArrayList<KnownGene>(this.knownGenes.get(currChrom));
				}
			
			Variant variant=new Variant();
			variant.ctx=ctx;
			variantStack.add(variant);
			
			
			
			//no more genes
			if(currListOfGenes.isEmpty())
				{
				while(!variantStack.isEmpty())
					{
					dump(w,variantStack.removeFirst());
					}
				continue;
				}
			
			//remove old variants
			while(  !variantStack.isEmpty() &&
					!currListOfGenes.isEmpty() &&
					variantStack.getFirst().ctx.getEnd()-1 < currListOfGenes.get(0).getTxStart()
					)
				{
				dump(w,variantStack.removeFirst());
				}
			
			//find genes to be processed
			int kgIndex=0;
			while(kgIndex< currListOfGenes.size())
				{
				KnownGene kg=currListOfGenes.get(kgIndex);
				if(kg.getTxEnd() < ctx.getStart())
					{
					challenge(variantStack,kg);
					currListOfGenes.remove(kgIndex);
					}
				else if(kg.getTxStart() > ctx.getEnd())
					{
					break;
					}
				else
					{
					kgIndex++;
					}
				}
			
			}
		while(!variantStack.isEmpty())
			{
			dump(w,variantStack.removeFirst());
			}
		}
	
	private void challenge(List<Variant> L,KnownGene gene) throws IOException
		{
		if(L.size()<2) return;
		List<Variant> variants=new ArrayList<Variant>(L);
		
		for(int i=1;i< variants.size();++i)
			{
			VariantContext ctx1=variants.get(i-1).ctx;
			VariantContext ctx2=variants.get(i).ctx;
			if(!ctx1.getChr().equals(ctx2.getChr())) throw new IllegalStateException();
			if(ctx1.getStart()> ctx2.getStart())
				{
				error("Data are not sorted");
				throw new IOException("VCF is not sorted");
				}
			}
		Map<Integer,Variant> genomicposzero2var=new HashMap<Integer,Variant>(variants.size());
		for(Variant var:variants)
			{
			if(genomicposzero2var.put(var.ctx.getStart()-1, var)!=null)
				{
				warning("duplicate variant at position "+gene.getChr()+":"+var.ctx.getStart());
				}
			}
		
			
		if(genomicSequence==null || !genomicSequence.getChrom().equals(gene.getChr()))
			{
			info("getting genomic Sequence for "+gene.getChr());
			genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, gene.getChr());
			}
				
			
		GeneticCode geneticCode=GeneticCode.getStandard();
    		
		StringBuilder wildRNA=new StringBuilder();
		MutedSequence mutRNA=new MutedSequence(wildRNA);
    					
		for(Variant var:variants)
			{
			var.position_in_cdna=-1;
			var.mutRNA=new MutedSequence(wildRNA);
			}
		if(gene.isPositiveStrand())
			{
	    	for(int exon_index=0;exon_index< gene.getExonCount();++exon_index)
	    		{
	    		KnownGene.Exon exon= gene.getExon(exon_index);
	    		for(int i= exon.getStart();
						i< exon.getEnd();
						++i)
					{
					if(i< gene.getCdsStart()) continue;
					if(i>=gene.getCdsEnd()) break;

					wildRNA.append(genomicSequence.charAt(i));
					
					Variant var=genomicposzero2var.get(i);
					if(var!=null)
						{
						var.position_in_cdna=wildRNA.length()-1;
						char mutBase=var.ctx.getAlternateAllele(0).getBaseString().charAt(0);
						mutRNA.put( wildRNA.length()-1, mutBase 	);
						var.mutRNA.put( wildRNA.length()-1, mutBase );
						}
						
					}
	    		}
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
					if(i>= gene.getCdsEnd()) continue;
					if(i<  gene.getCdsStart()) break;
	
				
					wildRNA.append(AcidNucleics.complement(genomicSequence.charAt(i)));
					
					
					Variant var=genomicposzero2var.get(i);
						
					if(var!=null)
						{
						var.position_in_cdna=wildRNA.length()-1;
						char mutBase=AcidNucleics.complement(var.ctx.getAlternateAllele(0).getBaseString().charAt(0));
						mutRNA.put( wildRNA.length()-1, mutBase 	);
						var.mutRNA.put( wildRNA.length()-1, mutBase );
						}
						
    				}
    			--exon_index;
    			}
			}
        	
    	for(Variant var:variants)
    		{
    		if(var.position_in_cdna==-1) continue;
    		int pos_aa=var.position_in_cdna/3;
    		int mod= var.position_in_cdna%3;
    		ProteinCharSequence wildProt=new ProteinCharSequence(geneticCode, wildRNA);
    		ProteinCharSequence mutProt=new ProteinCharSequence(geneticCode, mutRNA);
    		ProteinCharSequence varProt=new ProteinCharSequence(geneticCode, var.mutRNA);
    		char aa1=wildProt.charAt(pos_aa);
    		char aa2=mutProt.charAt(pos_aa);
    		char aa3=varProt.charAt(pos_aa);
    		if(aa3==aa2) continue;//already known mutation
    		if(aa3==aa1) continue;//silent mutation
    		int first_base_codon_in_cdna=var.position_in_cdna-mod;
    		StringBuilder sb=new StringBuilder();
    		sb.append(gene.getName());
    		sb.append("|");
    		sb.append(gene.isPositiveStrand()?'+':'-');
    		sb.append("|");
    		sb.append((1+var.position_in_cdna));
    		sb.append("|");
    		sb.append(wildRNA.substring(first_base_codon_in_cdna, first_base_codon_in_cdna+3));
    		sb.append("|");
    		sb.append(aa1);
    		sb.append("|");
    		sb.append(mutRNA.subSequence(first_base_codon_in_cdna, first_base_codon_in_cdna+3));
    		sb.append("|");
    		sb.append(aa2);
    		sb.append("|");
    		sb.append(var.mutRNA.subSequence(first_base_codon_in_cdna, first_base_codon_in_cdna+3));
    		sb.append("|");
    		sb.append(aa3);
    		sb.append("|");
    		for(Variant var2:variants)
    			{
    			if(var==var2) continue;
    			if(var2.position_in_cdna==-1) continue;
    			if(var2.position_in_cdna/3 != pos_aa) continue;
    			if(var2.ctx.hasID()) sb.append(var2.ctx.getID()+"_");
    			sb.append(var2.ctx.getStart());
    			sb.append("/");
    			}
    		
    		var.set.add(sb.toString());
    		//System.err.println(wildProt);
    		//System.err.println(mutProt);
    		}
			
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VCFStopCodon";
		}
	
	@Override
	public String getProgramDescription() {
		return  "Idea from @SolenaLS.";
		}
	
	@Override
	public void printOptions(PrintStream out) {
		out.println(" -R (file) indexed Fasta genome REFERENCE.");
		out.println(" -k (uri) KnownGene data URI/File. should look like"+ DEFAULT_KG_URI+"" +
				" . Beware chromosome names are formatted the same as your REFERENCE.");
		super.printOptions(out);
		}
	
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"R:k:"))!=-1)
			{
			switch(c)
				{
				case 'R': this.REF=new File(opt.getOptArg());break;
				case 'k': this.kgURI=opt.getOptArg();break;
				default:
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		if(this.REF==null)
			{
			error("Undefined REFERENCE.");
			return -1;
			}
		
		
		try
			{

			return super.doWork(opt.getOptInd(), args);
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
	
	public static void main(String[] args)
		{
		new VCFStopCodon().instanceMainWithExit(args);
		}
	}