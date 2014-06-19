package com.github.lindenb.jvarkit.tools.vcfdo;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Set;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.util.doid.DiseaseOntoglogyTree;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfFilterDoid
	extends AbstractVCFDiseaseOntology
	{
	
   /* list of DOID accessions for gene having a DOID-term children of the user output. */
	private Set<String> CHILD_OF=new HashSet<String>();
   
    
    /* Filter name */
    private String FILTER="DOID";
 
    
	    
	@Override
	public String getProgramDescription()
		{
		return "Set Filters for VCF annotated with SNPEFF or VEP tesing if Genes are part / are not part of a GeneOntology tree. ";
		}

	
	@Override
	protected void doWork(VcfIterator r, VariantContextWriter w)
			throws IOException
		{
		super.readDiseaseOntoglogyTree();
		Set<DiseaseOntoglogyTree.Term> positive_terms=new HashSet<DiseaseOntoglogyTree.Term>(CHILD_OF.size());
		for(String acn:CHILD_OF)
			{
			DiseaseOntoglogyTree.Term t=super.diseaseOntoglogyTree.getTermByAccession(acn);
			if(t==null)
				{
				throw new IOException("Cannot find GO acn "+acn);
				}
			positive_terms.add(t);
			}
		
		this.info("positive_terms:"+positive_terms);
	
		VCFHeader header=r.getHeader();

		
		super.readDiseaseOntoglogyAnnotations();
		super.loadEntrezGenes(header.getSequenceDictionary());
		

		
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine( new VCFFilterHeaderLine(
				FILTER
				,"Flag  DOID terms "+CHILD_OF+" using "+super.DOI_INPUT+" and "+super.DOI_ANN
				));
		
		
		w.writeHeader(h2);
	
		while(r.hasNext())
			{
			VariantContext ctx=r.next();
			VariantContextBuilder b=new VariantContextBuilder(ctx);

			
			Set<Integer> geneIds=getGeneIds(ctx);
		
			boolean set_positive=false;
			
			for(Integer GN:geneIds)
				{
				Set<DiseaseOntoglogyTree.Term> goset=super.gene2doid.get(GN);
				if(goset==null) continue;
				for(DiseaseOntoglogyTree.Term goParent:positive_terms)
					{
					for(DiseaseOntoglogyTree.Term t:goset)
						{
						if(goParent.hasDescendant(t.getAcn()))
							{
							set_positive=true;
							break;
							}
						}
					if(set_positive) break;
					}
				
				
				if(set_positive) break;
				}
			Set<String> filters=new HashSet<String>(ctx.getFilters());
			if(set_positive)
				{
				filters.add(FILTER);
				}
			
			b.filters(filters);

			w.add(b.make());
			}
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println("-C (term) ist of DOID accessions for gene having a DOID-term children of the user output");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "C:"))!=-1)
			{
			switch(c)
				{
				case 'C':CHILD_OF.add(opt.getOptArg()); break;
				default: 
					{
					switch(handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		System.err.println("I'm working on this. Sorry. Exiting");
		if("A".equals("A")) System.exit(-1);
		return doWork(opt.getOptInd(), args);
		}

	
	public static void main(String[] args)
		{
		new VcfFilterDoid().instanceMainWithExit(args);
		}
	}
