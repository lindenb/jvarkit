package com.github.lindenb.jvarkit.tools.vcfstripannot;

import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter;


public class VCFStripAnnotations extends AbstractVCFFilter
	{
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" Removes one or more field from the INFO column of a VCF.";
    
	@Option(shortName="K", doc="remove this INFO key",minElements=0)	
	public Set<String> KEY=new HashSet<String>();

	@Option(shortName="NF", doc="Reset the FILTER column",optional=true)	
	public boolean RESET_FILTER=false;

	@Override
	protected void doWork(VcfIterator r, VariantContextWriter w)
			throws IOException
			{
			VCFHeader header=r.getHeader();
			
			
			VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
			
			for(Iterator<VCFInfoHeaderLine> h=h2.getInfoHeaderLines().iterator();
					h.hasNext();)
				{
				VCFInfoHeaderLine vih=h.next();
				if(this.KEY.contains(vih.getID()))
					h.remove();
				}
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName(),"REMOVED: "+this.KEY.toString()));
			w.writeHeader(h2);
			
			while(r.hasNext())
				{
				VariantContext ctx=r.next();
				VariantContextBuilder b=new VariantContextBuilder(ctx);
				for(String key:KEY) b.rmAttribute(key);
				if(RESET_FILTER) b.unfiltered();
				w.add(b.make());
				}		
			}
	
	
	public static void main(String[] args) throws IOException
		{
		new VCFStripAnnotations().instanceMainWithExit(args);
		}
	}
