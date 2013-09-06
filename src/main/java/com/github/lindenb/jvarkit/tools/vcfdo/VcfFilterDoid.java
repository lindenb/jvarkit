package com.github.lindenb.jvarkit.tools.vcfdo;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.vcf.VcfIterator;

import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFFilterHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.util.doid.DiseaseOntoglogyTree;

public class VcfFilterDoid extends AbstractVCFDiseaseOntology
	{
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" Set Filters for VCF annotated with SNPEFF or VEP tesing if Genes are part / are not part of a GeneOntology tree. ";
	
    @Option(shortName="CHILD", doc="list of DOID accessions for gene having a DOID-term children of the user output.",minElements=0)
	public Set<String> CHILD_OF=new HashSet<String>();
   
    
    @Option(shortName="FILTER_NAME", doc="Filter name.",optional=true)
	public String FILTER="DOID";
 
    
    

	
	@Override
	public String getVersion()
		{
		return "1.0";
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
		
		LOG.info("positive_terms:"+positive_terms);
	
			
		
		super.readDiseaseOntoglogyAnnotations();
		super.loadEntrezGenes();
		
		VCFHeader header=r.getHeader();

		
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
	
	public static void main(String[] args)
		{
		new VcfFilterDoid().instanceMainWithExit(args);
		}
	}
