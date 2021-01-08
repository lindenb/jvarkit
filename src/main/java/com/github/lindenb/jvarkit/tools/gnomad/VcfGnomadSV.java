/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


*/
package com.github.lindenb.jvarkit.tools.gnomad;

import java.io.File;
import java.nio.file.Path;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.variant.sv.StructuralVariantComparator;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

# Example:

```
java -jar dist/vcfgnomadsv.jar \
	-g src/test/resources/gnomad_v2_sv.sites.vcf.gz \
	./src/test/resources/manta.B00GWGD.vcf.gz
```

END_DOC
 */
@Program(name="vcfgnomadsv",
	description="Peek annotations from gnomad structural variants",
	keywords={"vcf","annotation","gnomad","sv"},
	creationDate="20190814",
	modificationDate="20191002"
)
public class VcfGnomadSV extends Launcher{
	
	private static final Logger LOG = Logger.build(VcfGnomadSV.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-g","--gnomad"},description="Gnomad-SV VCF file. see https://gnomad.broadinstitute.org/downloads#structural-variants",required=true)
	private Path gnomadVcfSvPath = null;
	@Parameter(names={"-p","--prefix"},description="INFO field prefix")
	private String prefix = "GNOMAD_";
	@Parameter(names={"--in-gnomad-filter"},description="If not empty, set this FILTER is variant was found in gnomad")
	private String in_gnomad_filter = "";
	@Parameter(names={"--discordant_svtype"},description="If not empty, set this FILTER if SVTYPE are discordants")
	private String discordant_svtype_filter = "";
	@Parameter(names={"--any-overlap-filter"},description="If not empty, set this FILTER if any variant in gnomad is found overlaping the variant BUT we didn't find a correct match")
	private String any_overlap_filter = "";

	@ParametersDelegate
	private StructuralVariantComparator svComparator = new StructuralVariantComparator();

	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator iter,
			final VariantContextWriter out
			)
		{
		final VCFReader gnomadVcfReader = VCFReaderFactory.makeDefault().open(this.gnomadVcfSvPath,true);
		final VCFHeader gnomadHeader = gnomadVcfReader.getHeader();
		final SAMSequenceDictionary gnomadDict= SequenceDictionaryUtils.extractRequired(gnomadHeader);
		final ContigNameConverter contigNameConverter = ContigNameConverter.fromOneDictionary(gnomadDict);
		this.svComparator.setContigComparator((A,B)->{
			final String c1 = contigNameConverter.apply(A);
			if(StringUtils.isBlank(c1)) return false;
			final String c2 = contigNameConverter.apply(A);
			return c1.equals(c2);
		});
		
		final VCFHeader h0 = iter.getHeader();
		final SAMSequenceDictionary dict= h0.getSequenceDictionary();
		
		if(!SequenceDictionaryUtils.isGRCh37(dict)) {
			LOG.warn("Input is NOT GRCh37 ?");
			// can be lift over
			}
		
		

		final VCFHeader h2 = new VCFHeader(h0);
		if(!StringUtils.isBlank(this.in_gnomad_filter)) {
			h2.addMetaDataLine(new VCFFilterHeaderLine(this.in_gnomad_filter, "Variant is found in GNOMADSV"));
		}
		if(!StringUtils.isBlank(this.discordant_svtype_filter)) {
			h2.addMetaDataLine(new VCFFilterHeaderLine(this.discordant_svtype_filter, "SVTYPE are discordants."));
		}
		if(!StringUtils.isBlank(this.any_overlap_filter)) {
			h2.addMetaDataLine(new VCFFilterHeaderLine(this.any_overlap_filter, "A variant in gnomad is found overlaping the variant"));
		}

		for(final VCFInfoHeaderLine h: gnomadHeader.getInfoHeaderLines()) {
			final VCFInfoHeaderLine hd2 = VCFUtils.renameVCFInfoHeaderLine(h, this.prefix + h.getID());
			h2.addMetaDataLine(hd2);
		}
		for(final VCFFilterHeaderLine h: gnomadHeader.getFilterLines()) {
			final VCFFilterHeaderLine hd2 = new VCFFilterHeaderLine(this.prefix + h.getID(),h.getDescription());
			h2.addMetaDataLine(hd2);
		}
		
		final VCFInfoHeaderLine ghomadIdInfo = new VCFInfoHeaderLine(this.prefix + "ID",1,VCFHeaderLineType.String,"Original ID column in gnomad");
		h2.addMetaDataLine(ghomadIdInfo);
		
		JVarkitVersion.getInstance().addMetaData(this, h2);
		final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(dict).logger(LOG).build();
		
		out.writeHeader(h2);
		while(iter.hasNext()) {
			final VariantContext ctx = progress.apply(iter.next());
			
			final String ctg= contigNameConverter.apply(ctx.getContig());
			if(StringUtil.isBlank(ctg)) {
				out.add(ctx);
				continue;
				}
			
			
			VariantContext cgtx=null;
			final CloseableIterator<VariantContext> iter2 = gnomadVcfReader.query(
					ctg,
					Math.max(1,ctx.getStart()-this.svComparator.getBndDistance()),
					ctx.getEnd()+this.svComparator.getBndDistance()
					);
			boolean found_any_overlap = false;
			
			while(iter2.hasNext()) {
				final VariantContext ctx2 = iter2.next();
				
				if(this.svComparator.testSimpleOverlap(ctx, ctx2)) {
					found_any_overlap = true;
					}
				
				if(this.svComparator.test(ctx, ctx2)) {
					cgtx=ctx2;
					break;
					}
				}
			iter2.close();
			
			final VariantContextBuilder vcb=new VariantContextBuilder(ctx);

			
			if(cgtx==null) {
				if(found_any_overlap && !StringUtil.isBlank(this.any_overlap_filter)) {
					vcb.filter(this.any_overlap_filter);
					out.add(vcb.make());
					continue;
					}
				
				out.add(ctx);
				continue;
				}
			for(final String key: cgtx.getAttributes().keySet()) {
				vcb.attribute(this.prefix+key,cgtx.getAttribute(key));
				}
			for(final String f: cgtx.getFilters()) {
				vcb.filter(this.prefix+f);
				}
			
			if(cgtx.hasID()) {
				if(!ctx.hasID())
					{
					vcb.id(cgtx.getID());
					}
				vcb.attribute(ghomadIdInfo.getID(),cgtx.getID());
				}
			
			if(!StringUtil.isBlank(this.in_gnomad_filter)) {
				vcb.filter(this.in_gnomad_filter);
				}
			if(!StringUtil.isBlank(this.discordant_svtype_filter) && 
				!ctx.getAttributeAsString(VCFConstants.SVTYPE,".").equals(cgtx.getAttributeAsString(VCFConstants.SVTYPE, "."))) {
				vcb.filter(this.discordant_svtype_filter);
				}
			out.add(vcb.make());
			}
		out.close();
		progress.close();
		gnomadVcfReader.close();
		return 0;
		}

	public VcfGnomadSV() {
		}
	
	@Override
	public int doWork(final List<String> args) {
		return doVcfToVcf(args,this.outputFile);
		}

public static void main(final String[] args) {
	new VcfGnomadSV().instanceMainWithExit(args);
	}
}
