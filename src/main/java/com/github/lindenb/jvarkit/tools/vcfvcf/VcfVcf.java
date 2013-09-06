package com.github.lindenb.jvarkit.tools.vcfvcf;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.picard.vcf.VcfIterator;

import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter;



import org.broad.tribble.readers.TabixReader;


public class VcfVcf extends AbstractVCFFilter
	{
	 private static Log LOG=Log.getInstance(VcfVcf.class); 
	
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"Get the INFO from a VCF and use it for another VCF. ";

	@Option(shortName="TBX",doc="The VCF file indexed with TABIX. Source of the annotations")
	public String TABIX;
	
	@Option(shortName="INFO",doc="The INFO keys to grab.",minElements=0)	
	public Set<String> INFO_IDS=new LinkedHashSet<String>();
	
	
	
	
	private TabixReader tabixReader =null;
	
	@Option(shortName="RIF",doc="Replace the INFO field if it exists.",minElements=0)	
	public boolean REPLACE_INFO_FIELD=true;
	@Option(shortName="RID",doc="Replace the ID field if it exists.",optional=true)			
	public boolean REPLACE_ID=true;
	@Option(shortName="RAM",doc="REF allele matters.",optional=true)			
	public boolean REF_ALLELE_MATTERS=true;
	@Option(shortName="AAM",doc="ALT alleles matters.",optional=true)			
	public boolean ALT_ALLELES_MATTERS=false;
	@Option(shortName="ACF",doc="Flag to set if alternate alleles conflict.",optional=true)			
	public String ALT_CONFLICt_FLAG=null;
	
	
	@Override
	protected void doWork(VcfIterator r, VariantContextWriter w)
			throws IOException
		{
		VCFCodec codeIn3=new VCFCodec();	
		String line;
		
		StringWriter sw=new StringWriter();
		LOG.info("opening "+this.TABIX);
	    TabixReader tabix= new TabixReader(this.TABIX);
		while((line=tabix.readLine())!=null)
			{
			if(!line.startsWith(VCFHeader.HEADER_INDICATOR))
				{
				break;
				}
			sw.append(line).append("\n");
			}
		VCFHeader header3=new VcfIterator(new ByteArrayInputStream(sw.toString().getBytes())).getHeader();
		VCFHeader header1=r.getHeader();
		
		VCFHeader h2=new VCFHeader(header1.getMetaDataInInputOrder(),header1.getSampleNamesInOrder());
		for(String infoId:this.INFO_IDS)
			{
			VCFInfoHeaderLine vihl=header3.getInfoHeaderLine(infoId);
			if(vihl==null)
				{
				LOG.warn("Not INFO="+infoId+" in "+TABIX);
				continue;
				}
			if(h2.getInfoHeaderLine(infoId)!=null)
				{
				LOG.warn("Input already contains INFO="+vihl);
				}
			h2.addMetaDataLine(vihl);
			}
		
		w.writeHeader(h2);
		while(r.hasNext())
			{
			VariantContext ctx1=r.next();
			VariantContextBuilder  vcb=new VariantContextBuilder(ctx1);
			String line2;
			String BEST_ID=null;
			boolean best_id_match_alt=false;

			List<VariantContext> variantsList=new ArrayList<VariantContext>();
			
			TabixReader.Iterator iter=tabixReader.query(ctx1.getChr()+":"+ctx1.getStart()+"-"+ctx1.getEnd());
			while(iter!=null && (line2=iter.next())!=null)
				{
				VariantContext ctx3=codeIn3.decode(line2);
				if(ctx3.getStart()!=ctx1.getStart()) continue;
				if(ctx3.getEnd()!=ctx1.getEnd()) continue;
				
				if( ctx1.getReference().equals(ctx3.getReference()) &&
					ctx1.getAlternateAlleles().equals(ctx3.getAlternateAlleles())
					)
					{
					variantsList.clear();
					variantsList.add(ctx3);
					break;
					}
				else
					{
					variantsList.add(ctx3);
					}
				}
			for(VariantContext ctx3:variantsList)
				{
				if(ctx3.getID()!=null && this.REPLACE_ID)
					{
					if(BEST_ID!=null && best_id_match_alt)
						{
						continue;
						}
					BEST_ID=ctx3.getID();
					best_id_match_alt=ctx1.getAlternateAlleles().equals(ctx3.getAlternateAlleles());
					}
				
				if(this.REF_ALLELE_MATTERS && !ctx1.getReference().equals(ctx3.getReference()))
					{
					continue;
					}
				if(this.ALT_ALLELES_MATTERS && !ctx1.getAlternateAlleles().equals(ctx3.getAlternateAlleles()))
					{
					continue;
					}
				for(String id:this.INFO_IDS)
					{
					Object info3=ctx1.getAttribute(id);
					if(info3==null) continue;
					
					Object info1=ctx1.getAttribute(id);
					if(info1!=null && !this.REPLACE_INFO_FIELD)
						{
						continue;
						}
					vcb.attribute(id, info3);
					}
				
				if(ALT_CONFLICt_FLAG!=null && !ctx1.getAlternateAlleles().equals(ctx3.getAlternateAlleles()))
					{
					vcb.attribute(ALT_CONFLICt_FLAG, true);
					}
				
				}
			if(BEST_ID!=null)
				{
				vcb.id(BEST_ID);
				}
			w.add(vcb.make());
			}
		tabix.close();
		}
	
	
	public static void main(String[] args) throws IOException
		{
		new VcfVcf().instanceMainWithExit(args);
		}
}
