package com.github.lindenb.jvarkit.tools.vcfgo;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import htsjdk.samtools.cmdline.Option;
import htsjdk.samtools.cmdline.Usage;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.util.go.GoTree;

public class VcfFilterGo extends AbstractVcfGeneOntology
	{
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" Set Filters for VCF annotated with SNPEFF or VEP tesing if Genes are part / are not part of a GeneOntology tree. ";
	
    @Option(shortName="CHILD", doc="list of GO accessions for gene having a GO-term children of the user output.",minElements=0)
	public Set<String> CHILD_OF=new HashSet<String>();
   
    
    @Option(shortName="FILTER_NAME", doc="Filter name.",optional=true)
	public String FILTER="GO";
 
    
    

	
	@Override
	public String getVersion()
		{
		return "1.0";
		}

	@Override
	protected void doWork(VcfIterator r, VariantContextWriter w)
			throws IOException
		{
		super.readGO();
		Set<GoTree.Term> positive_terms=new HashSet<GoTree.Term>(CHILD_OF.size());
		for(String acn:CHILD_OF)
			{
			GoTree.Term t=super.goTree.getTermByAccession(acn);
			if(t==null)
				{
				throw new IOException("Cannot find GO acn "+acn);
				}
			positive_terms.add(t);
			}
		
		LOG.info("positive_terms:"+positive_terms);
	
		
		VCFHeader header=r.getHeader();
		
		
		
		super.readGOA();
		super.loadBiomartHGNC();
		
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine( new VCFFilterHeaderLine(
				FILTER
				,"Flag  GO terms "+CHILD_OF+" using "+super.GO+" and "+super.GOA
				));
		
		
		w.writeHeader(h2);
	
		while(r.hasNext())
			{
			VariantContext ctx=r.next();
			VariantContextBuilder b=new VariantContextBuilder(ctx);

			
			Set<String> geneNames=getGeneNames(ctx);
		
			boolean set_positive=false;
			
			for(String GN:geneNames)
				{
				Set<GoTree.Term> goset=super.name2go.get(GN);
				if(goset==null) continue;
				for(GoTree.Term goParent:positive_terms)
					{
					for(GoTree.Term t:goset)
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
	
	public static void main(String[] args)
		{
		new VcfFilterGo().instanceMainWithExit(args);
		}
	}
