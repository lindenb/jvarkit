/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

END_DOC
 */
@Program(name="vcfgnomadsv",
	description="Peek annotations from gnomad structural variants",
	keywords={"vcf","annotation","gnomad","sv"},
	creationDate="20190814",
	modificationDate="20190814"
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
	@Parameter(names={"--fraction"},description="two segments are identical if they overlap at fraction 'x' ")
	private double fraction_overlap = 0.7;
	@Parameter(names={"--bnd-distance"},description="two bnd are identical if they're distant from  'x' bases")
	private int bnd_distance=10;



	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator iter,
			final VariantContextWriter out
			)
		{
		final VCFFileReader gnomadVcfReader = new VCFFileReader(this.gnomadVcfSvPath,true);
		final VCFHeader gnomadHeader = gnomadVcfReader.getFileHeader();
		final SAMSequenceDictionary gnomadDict= SequenceDictionaryUtils.extractRequired(gnomadHeader);
		final ContigNameConverter contigNameConverter = ContigNameConverter.fromOneDictionary(gnomadDict);
		
		
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
			final String svType = ctx.getAttributeAsString(VCFConstants.SVTYPE,"");
			if(StringUtil.isBlank(svType)) {
				LOG.warn("No SV Type in "+ctx.getContig()+":"+ctx.getStart());
				out.add(ctx);
				continue;
				}
			
			VariantContext cgtx=null;
			final CloseableIterator<VariantContext> iter2 = gnomadVcfReader.query(
					ctg,
					Math.max(1,ctx.getStart()-bnd_distance),
					ctx.getEnd()+bnd_distance
					);
			boolean found_any_overlap = false;
			
			while(iter2.hasNext()) {
				final VariantContext ctx2 = iter2.next();
				final String gSvType = ctx2.getAttributeAsString(VCFConstants.SVTYPE,"");
				if(StringUtil.isBlank(gSvType)) {
					LOG.warn("No SV Type in gnomad: "+ctx2.getContig()+":"+ctx2.getStart());
					continue;
					}
				if(CoordMath.overlaps(ctx.getStart(), ctx.getEnd(), ctx2.getStart(), ctx2.getEnd())) /* don't use 'ctx.overlap(ctx2)' chr1 vs 1 */{
					found_any_overlap = true;
				}
				
				if(svType.equals("BND") && gSvType.equals(svType))
					{
					if(Math.abs(ctx.getStart()-ctx2.getStart()) > this.bnd_distance) continue;
					if(Math.abs(ctx.getEnd()-ctx2.getEnd()) > this.bnd_distance) continue;
					cgtx=ctx2;
					break;
					}
				else if(svType.equals("BND") || gSvType.equals("BND")) {
					//ignore
					}
				else
					{
					if(!CoordMath.overlaps(ctx.getStart(), ctx.getEnd(), ctx2.getStart(), ctx2.getEnd())) /* don't use 'ctx.overlap(ctx2)' chr1 vs 1 */{
						continue;
						}
					final int left= Math.max(ctx.getStart(),ctx2.getStart());
					final int right= Math.min(ctx.getEnd(),ctx2.getEnd());
					if(left > right) {
						LOG.debug("strange overlap??\n"+
								ctx+"\n"+
								ctx2);
						continue;
						}
					double len = CoordMath.getLength(left, right);
					//System.err.println("#");
					//System.err.println(len);
					//System.err.println(len/ctx.getLengthOnReference());
					//System.err.println(len/ctx2.getLengthOnReference());
					
					if(len/ctx.getLengthOnReference() < this.fraction_overlap ) continue;
					if(len/ctx2.getLengthOnReference() < this.fraction_overlap ) continue;
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
			if(!ctx.hasID())
				{
				vcb.id(cgtx.getID());
				}
			if(!StringUtil.isBlank(this.in_gnomad_filter)) {
				vcb.filter(this.in_gnomad_filter);
				}
			if(!StringUtil.isBlank(this.discordant_svtype_filter) && !svType.equals(cgtx.getAttributeAsString(VCFConstants.SVTYPE, ""))) {
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
