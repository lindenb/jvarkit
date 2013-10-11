package com.github.lindenb.jvarkit.util.vcf;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.CommonInfo;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypesContext;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFContigHeaderLine;
import org.broadinstitute.variant.vcf.VCFFilterHeaderLine;
import org.broadinstitute.variant.vcf.VCFFormatHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;


public class XMLVcfWriter implements VariantContextWriter
	{
	private XMLStreamWriter writer;
	
	
	public XMLVcfWriter(XMLStreamWriter writer)
		{
		this.writer=writer;
		}		
	protected void start(String tag) throws XMLStreamException
		{
		this.writer.writeStartElement(tag);
		}
	
	protected void attribute(String tag,Object o) throws XMLStreamException
		{
		this.writer.writeAttribute(tag, String.valueOf(o));
		}
	
	protected void end() throws XMLStreamException
		{
		this.writer.writeEndElement();
		}
	
	protected void element(String tag,Object content) throws XMLStreamException
		{
		if(content==null)
			{
			this.writer.writeEmptyElement(tag);
			return;
			}
		start(tag);
		characters(content);
		end();
		}

	protected void characters(Object content)
			throws XMLStreamException
		{
		if(content==null) return;
		this.writer.writeCharacters(String.valueOf(content));
		}
	
	@Override
	public void writeHeader(VCFHeader header)
		{
		try
			{
			start("vcf");
			start("header");
			start("infos");
			for(VCFInfoHeaderLine h:header.getInfoHeaderLines())
				{
				start("info");
				attribute("ID", h.getID());
				element("key",h.getKey());
				element("count",h.getCount());
				element("description",h.getDescription());
				element("countType",h.getCountType());
				element("value",h.getValue());
				h.getCountType();
				end();
				}
			end();
			
			start("formats");
			for(VCFFormatHeaderLine h: header.getFormatHeaderLines())
				{
				start("format");
				attribute("ID", h.getID());
				element("key",h.getKey());
				element("count",h.getCount());
				element("description",h.getDescription());
				element("countType",h.getCountType());
				element("value",h.getValue());
				end();
				}
			end();
			
			start("filters");
			for(VCFFilterHeaderLine h: header.getFilterLines())
				{
				start("filter");
				attribute("ID", h.getID());
				element("key",h.getKey());
				element("value",h.getValue());
				end();
				}
			end();
			
			
			start("contigs");
			for(VCFContigHeaderLine h:header.getContigLines())
				{
				start("contig");
				attribute("ID", h.getID());
				attribute("tid", h.getContigIndex());
				element("key",h.getKey());
				element("value",h.getValue());
				end();
				}
			end();
			
			start("samples");
			for(String name:header.getSampleNamesInOrder())
				{
				start("sample");
				attribute("index", header.getSampleNameToOffset().get(name));
				characters(name);
				end();
				}
			end();
			
			for(VCFHeaderLine meta:header.getMetaDataInInputOrder())
				{
				
				}
			
			end();//header
			start("variations");
			}
		catch (XMLStreamException e)
			{
			// TODO: handle exception
			}
		}
	
    // should we write genotypes or just sites?
	private boolean doNotWriteGenotypes=false;
	
	@Override
	public void add(VariantContext variant)
		{
		try
			{
			 if ( doNotWriteGenotypes )
				 variant = new VariantContextBuilder(variant).noGenotypes().make();
			
			 Map<Allele, String> alleleMap = buildAlleleMap(variant);
			 
			start("variation");
			element("chrom",variant.getChr());
			element("start",variant.getStart());
			element("end",variant.getEnd());
			element("id",variant.getID());
			element("ref",variant.getReference().getDisplayString());
			
			if ( variant.isVariant() )
				{
				for(Allele a:variant.getAlternateAlleles())
					{
					element("alt",a.getDisplayString());
					}
				}
			
			if(variant.hasLog10PError())
				{
				element("qual", variant.getLog10PError());
				}
			
			
			start("filters");
			if(variant.isFiltered())
				{
				for(String s: variant.getFilters())
					{
					element("filter",s);
					}
				}
			else if(variant.filtersWereApplied())
				{
				element("filter",VCFConstants.PASSES_FILTERS_v4);
				}
			end();
				
			start("infos");
			CommonInfo info=variant.getCommonInfo();
			for(String key:info.getAttributes().keySet())
				{
				start("info");
				attribute("key", key);
				Object v=info.getAttribute(key);
				
				end();
				}
			end();
			
			final GenotypesContext gc = variant.getGenotypes();
			
			
			if(variant.hasGenotypes())
				{
				start("genotypes");
				for(String sample:variant.getSampleNames())
					{
					Genotype g=variant.getGenotype(sample);
					if(g==null) continue;
					start("genotype");
					attribute("available",g.isAvailable());
					attribute("called",g.isCalled());
					attribute("het",g.isHet());
					attribute("hom",g.isHom());
					attribute("homRef",g.isHomRef());
					attribute("homVar",g.isHomVar());
					attribute("mixed",g.isMixed());
					attribute("noCall",g.isNoCall());
					attribute("nonInformative",g.isNonInformative());
					attribute("filtered",g.isFiltered());
					attribute("phased",g.isPhased());
					element("sample",g.getSampleName());
					if(g.hasAD())
						{
						start("ad");
						for(int ad:g.getAD())
							{
							element("value", ad);
							}
						g.getAD();
						end();
						}
					if(g.hasDP())
						{
						element("dp", g.getDP());
						}
					if(g.hasGQ())
						{
						element("gq", g.getGQ());
						}
					if(g.hasPL())
						{
						start("pl");
						for(int v:g.getPL())
							{
							element("value", v);
							}
						g.getAD();
						end();
						}
					if(g.hasLikelihoods())
						{
						
						}
					
					start("alleles");
					for(Allele a:g.getAlleles())
						{
						element("allele",a.getBaseString());
						}
					end();
					
					
					end();
					}
				end();
				}
			end();
			}
		catch (XMLStreamException e)
			{
			e.printStackTrace();
			}
		}
	List<String> x(VariantContext vc,VCFHeader header)
		{
		Set<String> keys = new HashSet<String>();

		 boolean sawGoodGT = false;
	        boolean sawGoodQual = false;
	        boolean sawGenotypeFilter = false;
	        boolean sawDP = false;
	        boolean sawAD = false;
	        boolean sawPL = false;
	        for ( final Genotype g : vc.getGenotypes() ) {
	            keys.addAll(g.getExtendedAttributes().keySet());
	            if ( g.isAvailable() ) sawGoodGT = true;
	            if ( g.hasGQ() ) sawGoodQual = true;
	            if ( g.hasDP() ) sawDP = true;
	            if ( g.hasAD() ) sawAD = true;
	            if ( g.hasPL() ) sawPL = true;
	            if (g.isFiltered()) sawGenotypeFilter = true;
	        }

	        if ( sawGoodQual ) keys.add(VCFConstants.GENOTYPE_QUALITY_KEY);
	        if ( sawDP ) keys.add(VCFConstants.DEPTH_KEY);
	        if ( sawAD ) keys.add(VCFConstants.GENOTYPE_ALLELE_DEPTHS);
	        if ( sawPL ) keys.add(VCFConstants.GENOTYPE_PL_KEY);
	        if ( sawGenotypeFilter ) keys.add(VCFConstants.GENOTYPE_FILTER_KEY);

	        List<String> sortedList = ParsingUtils.sortList(new ArrayList<String>(keys));

	        // make sure the GT is first
	        if ( sawGoodGT ) {
	            List<String> newList = new ArrayList<String>(sortedList.size()+1);
	            newList.add(VCFConstants.GENOTYPE_KEY);
	            newList.addAll(sortedList);
	            sortedList = newList;
	        }

	        if ( sortedList.isEmpty() && header.hasGenotypingData() ) {
	            // this needs to be done in case all samples are no-calls
	            return Collections.singletonList(VCFConstants.GENOTYPE_KEY);
	        } else {
	            return sortedList;
	        }
		}
	
	@Override
	public void close()
		{
		try
			{
			end();//variations
			end();//vcf
			writer.close();
			}
		catch (XMLStreamException e)
			{
			e.printStackTrace();
			}
		
		}
	
    private static Map<Allele, String> buildAlleleMap(final VariantContext vc) {
    final Map<Allele, String> alleleMap = new HashMap<Allele, String>(vc.getAlleles().size()+1);
    alleleMap.put(Allele.NO_CALL, VCFConstants.EMPTY_ALLELE); // convenience for lookup

    final List<Allele> alleles = vc.getAlleles();
    for ( int i = 0; i < alleles.size(); i++ ) {
        alleleMap.put(alleles.get(i), String.valueOf(i));
    }

    return alleleMap;
}

	
	}
