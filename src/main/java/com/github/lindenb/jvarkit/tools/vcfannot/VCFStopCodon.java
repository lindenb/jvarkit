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
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Pattern;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeaderLineCount;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamSequenceRecordTreeMap;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;





/**
 * 
 *
 */
public class VCFStopCodon extends AbstractVCFFilter2
	{
	private static final String DEFAULT_KG_URI="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz";
	private File REF=null;
	private String kgURI=DEFAULT_KG_URI;

	private SamSequenceRecordTreeMap<KnownGene> knownGenes=null;
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
		int n_genes=0;
		this.knownGenes=new SamSequenceRecordTreeMap<KnownGene>(
				this.indexedFastaSequenceFile.getSequenceDictionary()
				);
		info("loading genes");
		Set<String> unknown=new HashSet<String>();
		BufferedReader in=IOUtils.openURIForBufferedReading(this.kgURI);
		String line;
		Pattern tab=Pattern.compile("[\t]");
		while((line=in.readLine())!=null)
			{
			if(line.isEmpty()) continue;
			String tokens[]=tab.split(line);
			KnownGene g=new KnownGene(tokens);

			if(!this.knownGenes.put(g.getChr(), g.getTxStart()+1, g.getTxEnd(),g))
				{
				if(!unknown.contains(g.getChr()))
					{
					warning("The reference "+REF+" doesn't contain chromosome "+g.getChr());
					unknown.add(g.getChr());
					}
				}
			else
				{
				++n_genes;
				}
			}
		in.close();
		info("genes:"+n_genes);
		}
	
	
	
	private static class Variant
		{
		VariantContext ctx;
		int position_in_cdna=-1;
		MutedSequence mutRNA=null;
		Set<String> set=new HashSet<String>();
		}
	static final String TAG="STOPCODON";
	
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

		
        LinkedList<Variant> variantStack=new LinkedList<Variant>();
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
		while(r.hasNext())
			{
			VariantContext ctx=r.next();
			progress.watch(ctx.getChr(), ctx.getStart());
			if(!ctx.isSNP() || ctx.getAlternateAlleles().size()!=1)
				{
				w.add(ctx);
				continue;
				}
			Variant variant=new Variant();
			variant.ctx=ctx;
			variantStack.add(variant);
			
			if(variantStack.size()>3)
				{
				dump(w,variantStack.removeFirst());
				}
			
			if(variantStack.size()>2)
				{
				int prev_idx=variantStack.size()-3;
				if(!(variantStack.get(prev_idx).ctx.getChr().equals(ctx.getChr()) &&
					 variantStack.get(prev_idx).ctx.getStart()+2!=ctx.getStart()))
					{
					dump(w,variantStack.remove(prev_idx));
					}
				}
			
			if(variantStack.size()>1)
				{
				int prev_idx=variantStack.size()-2;
				if(!(variantStack.get(prev_idx).ctx.getChr().equals(ctx.getChr()) &&
					 variantStack.get(prev_idx).ctx.getStart()+1!=ctx.getStart()))
					{
					dump(w,variantStack.remove(prev_idx));
					}
				}
			if(variantStack.size()>1)
				{
				check(variantStack);
				}
			}
		while(!variantStack.isEmpty())
			{
			dump(w,variantStack.removeFirst());
			}
		}
	
	private void check(LinkedList<Variant> variants)
		{
		
		String chrom=variants.getFirst().ctx.getChr();
		int start1=variants.getFirst().ctx.getStart();
		int end1=variants.getLast().ctx.getStart();
		info(chrom+":"+start1+"-"+end1);
		
		List<KnownGene> genes=this.knownGenes.getOverlapping(
				chrom, start1,end1 //1-based
				);
		if(genes==null || genes.isEmpty())
			{
			return;
			}
			
		if(genomicSequence==null || !genomicSequence.getChrom().equals(chrom))
			{
			info("getting genomic Sequence for "+chrom);
			genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, chrom);
			}
				
		for(KnownGene gene:genes)
			{
			if(gene.isNonCoding()) continue;
			GeneticCode geneticCode=GeneticCode.getStandard();
    		
    		
	
    		StringBuilder wildRNA=new StringBuilder();
    		MutedSequence mutRNA=new MutedSequence(wildRNA);
    					

			for(Variant var:variants)
				{
				var.position_in_cdna=-1;
				var.mutRNA=new MutedSequence(wildRNA);
				}
			
        	for(int exon_index=0;exon_index< gene.getExonCount();++exon_index)
        		{
        		KnownGene.Exon exon= gene.getExon(exon_index);
        		for(int i= exon.getStart();
    					i< exon.getEnd();
    					++i)
    				{
    				if(i< gene.getCdsStart()) continue;
    				if(i>=gene.getCdsEnd()) break;
    				
    				for(Variant var:variants)
    					{
    					if(var.ctx.getStart()==i)
    						{
    						var.position_in_cdna=wildRNA.length();
    						}
    					}
    				
    				if(gene.isPositiveStrand())
    					{
    					wildRNA.append(genomicSequence.charAt(i));
    					}
    				else
    					{
    					wildRNA.insert(0,AcidNucleics.complement(genomicSequence.charAt(i)));
    					}
    				
    				for(Variant var:variants)
    					{
    					if(var.ctx.getStart()==i)
    						{
    						char mutBase=var.ctx.getAlternateAllele(0).getBaseString().charAt(0);
    						if(gene.isNegativeStrand())
    							{
    							mutBase=AcidNucleics.complement(mutBase);
    							}
    						mutRNA.put( wildRNA.length()-1, mutBase 	);
							var.mutRNA.put( wildRNA.length()-1, mutBase );
    						}
    					}
    				}
        		}
        	
        	for(Variant var:variants)
        		{
        		if(var.position_in_cdna==-1) continue;
        		int pos_aa=var.position_in_cdna/3;
        		//int mod= var.position_in_cdna%3;
        		ProteinCharSequence wildProt=new ProteinCharSequence(geneticCode, wildRNA);
        		ProteinCharSequence mutProt=new ProteinCharSequence(geneticCode, mutRNA);
        		ProteinCharSequence varProt=new ProteinCharSequence(geneticCode, var.mutRNA);
        		char aa1=wildProt.charAt(pos_aa);
        		char aa2=mutProt.charAt(pos_aa);
        		char aa3=varProt.charAt(pos_aa);
        		if(aa3==aa2) continue;//already known mutation
        		if(aa3==aa1) continue;//silent mutation
        		var.set.add("FOUUUNNNNNNNNNNNNNNNNNNNNNNNNNNNND");
        		}            	
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