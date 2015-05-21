package com.github.lindenb.jvarkit.tools.vcfvcf;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import com.github.lindenb.jvarkit.util.picard.cmdline.Option;
import com.github.lindenb.jvarkit.util.picard.cmdline.Usage;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter;
import com.github.lindenb.jvarkit.util.vcf.IndexedVcfFileReader;

@Deprecated
public class VcfVcf extends AbstractVCFFilter
	{
	 private static Log LOG=Log.getInstance(VcfVcf.class); 
	
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"Get the INFO from a VCF and use it for another VCF. ";

	@Option(shortName="TBX",doc="The VCF file indexed with TABIX. Source of the annotations")
	public String TABIX;
	
	@Option(shortName="INFO",doc="The INFO keys to grab.",minElements=0)	
	public Set<String> INFO_IDS=new LinkedHashSet<String>();
	
	
	
	
	@Option(shortName="RIF",doc="Replace the INFO field if it exists.",minElements=0)	
	public boolean REPLACE_INFO_FIELD=true;
	@Option(shortName="RID",doc="Replace the ID field if it exists.",optional=true)			
	public boolean REPLACE_ID=true;
	@Option(shortName="RAM",doc="REF allele matters.",optional=true)			
	public boolean REF_ALLELE_MATTERS=true;
	@Option(shortName="AAM",doc="ALT alleles matters.",optional=true)			
	public boolean ALT_ALLELES_MATTERS=false;
	@Option(shortName="ACF",doc="Flag to set if alternate alleles conflict.",optional=true)			
	public String ALT_CONFLICT_FLAG=null;
	
	
	@Override
	protected void doWork(VcfIterator r, VariantContextWriter w)
			throws IOException
		{
		CloseableIterator<VariantContext> iter=null;
		LOG.info("opening file: "+this.TABIX);
	    IndexedVcfFileReader tabix= new IndexedVcfFileReader(this.TABIX);
		VCFHeader header3=tabix.getHeader();
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
		
		if(ALT_CONFLICT_FLAG!=null)
			{
			h2.addMetaDataLine(new VCFInfoHeaderLine(ALT_CONFLICT_FLAG,1,VCFHeaderLineType.Flag,"conflict ALT allele with "+this.TABIX));
			}
		
		w.writeHeader(h2);
		while(r.hasNext())
			{
			VariantContext ctx1=r.next();
			
			VariantContextBuilder  vcb=new VariantContextBuilder(ctx1);
			String BEST_ID=null;
			boolean best_id_match_alt=false;

			List<VariantContext> variantsList=new ArrayList<VariantContext>();
			
		
			iter=tabix.iterator(ctx1.getChr(),
					Math.max(0,ctx1.getStart()-1),
					(ctx1.getEnd()+1)
					);
			
			while(iter.hasNext())
				{
				VariantContext ctx3=iter.next();
				if(!ctx3.getChr().equals(ctx1.getStart())) continue;
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
			CloserUtil.close(iter);iter=null;
			
			for(VariantContext ctx3:variantsList)
				{
				
				
				if(this.REF_ALLELE_MATTERS && !ctx1.getReference().equals(ctx3.getReference()))
					{
					continue;
					}
				if(this.ALT_ALLELES_MATTERS && !ctx1.getAlternateAlleles().equals(ctx3.getAlternateAlleles()))
					{
					continue;
					}
				
				if(ctx3.getID()!=null && this.REPLACE_ID)
					{
					if(BEST_ID!=null && best_id_match_alt)
						{
						//nothing
						}
					else
						{
						BEST_ID=ctx3.getID();
						best_id_match_alt=ctx1.getAlternateAlleles().equals(ctx3.getAlternateAlleles());
						}
					}
				
				
				for(String id:this.INFO_IDS)
					{
					Object info3=ctx3.getAttribute(id);
					if(info3==null)
						{
						continue;
						}
					Object info1=ctx1.getAttribute(id);
					if(info1!=null && !this.REPLACE_INFO_FIELD)
						{
						continue;
						}
					
					vcb.attribute(id, info3);
					}
				
				if(ALT_CONFLICT_FLAG!=null && !ctx1.getAlternateAlleles().equals(ctx3.getAlternateAlleles()))
					{
					vcb.attribute(ALT_CONFLICT_FLAG, true);
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
