/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.Closeable;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.variant.sv.OverlapComparator;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFReader;

/**

BEGIN_DOC

* For Structural Variants (identified with INFO/SVTYPE) we just compare the fraction of common overlap.
* For SNVs : spanning deletions (`*`) are not considered when comparing alleles

See also:

 * bcftools isec

END_DOC
*/

@Program(name="vcfin",
	description="Only prints variants that are contained/not contained into another VCF",
	keywords={"vcf","compare"},
	modificationDate="20250903",
	creationDate="20140204",
	biostars={287815}
	)
public class VcfIn extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.of(VcfIn.class);
	private enum Matcher {id,overlap,chrom_pos,chrom_pos_ref,chrom_pos_ref_any_alt,chrom_pos_ref_all_alt}
	
	@Parameter(names={"-i","--inverse"},description="Print variants that are not part of the VCF-database.")
	private boolean inverse = false;
	@Parameter(names={"--vcf","--vcfs","--file-list"},description= "External indexed VCFs. a file with the suffix '.list' is interpreted as a file containing the path to the indexed VCFs")
	private List<String> databaseVCFPaths = new ArrayList<>();
	@Parameter(names={"--filter"},description="FILTER name for SOFT filtering")
	private String filterName = "";
	@Parameter(names={"--sv-fraction"},description="How to match two SV. " + OverlapComparator.OPT_DESC,splitter = NoSplitter.class,converter = OverlapComparator.StringConverter.class)
	private OverlapComparator overlap_comparator = OverlapComparator.makeDefault();
	@Parameter(names={"--matcher"},description="How to match two SNVs. Note that Spanning deletion are ignored.")
	private Matcher matcher = Matcher.chrom_pos_ref_any_alt;

	
	private final List<DatabaseVcf> databaseVCFS = new ArrayList<>();
	
	
	private static class DatabaseVcf
	implements Closeable
		{
		final Path file;
		final VCFReader vcfReader;
		final ContigNameConverter contigNameConverter;
		DatabaseVcf(Path file) {
			this.file = file;
			this.vcfReader =VCFReaderFactory.makeDefault().open(this.file,true);
			final VCFHeader header = this.vcfReader.getHeader();
			final SAMSequenceDictionary dict = header.getSequenceDictionary();
			this.contigNameConverter = dict==null
					? ContigNameConverter.getIdentity()
					:ContigNameConverter.fromOneDictionary(dict);
			}
		
		
		@Override
		public void close() {
			try {vcfReader.close();}
			catch(IOException err) {LOG.error(err);}
			}
		}

	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	/** if dict are different , ignore contig*/
	private static boolean compare(
			final Matcher matcher,
			final VariantContext vc1,
			final VariantContext vc2
			)
		{
		switch(matcher) {
			case overlap:
				return true;//assume we already know they overlap
			case id:
				{
				if(!vc1.hasID() || !vc2.hasID()) return false;
				final Set<String> id1s =CharSplitter.of(VCFConstants.ID_FIELD_SEPARATOR).splitAsSet((vc1.getID()));
				final Set<String> id2s =CharSplitter.of(VCFConstants.ID_FIELD_SEPARATOR).splitAsSet((vc2.getID()));
				for(String id: id1s) {
					if(id.equals(VCFConstants.EMPTY_ID_FIELD)) continue;
					if(id2s.contains(id)) return true;
					}
				return false;
				}
			case chrom_pos:
				return vc1.getStart() == vc2.getStart();
			case chrom_pos_ref:
				{
				if(vc1.getStart() != vc2.getStart()) return false;
				if(vc1.getNAlleles()<1 || vc2.getNAlleles()<1) return false;
				return vc1.getReference().equals(vc2.getReference());
				}
			case chrom_pos_ref_any_alt:
				{
				if(!compare(Matcher.chrom_pos_ref,vc1,vc2)) return false;
				for(Allele alt: vc1.getAlternateAlleles()) {
					if(alt.equals(Allele.SPAN_DEL)) continue;
					if(vc2.hasAllele(alt)) return true;
					}
				return false;
				}
			case chrom_pos_ref_all_alt:
				{
				if(!compare(Matcher.chrom_pos_ref,vc1,vc2)) return false;
				if(vc1.getNAlleles()<2) return false;
				if(vc2.getNAlleles()<2) return false;
				boolean ok=false;
				for(Allele alt: vc1.getAlternateAlleles()) {
					if(alt.equals(Allele.SPAN_DEL)) continue;
					if(!vc2.hasAllele(alt)) return false;
					ok=true;
					}
				return ok;
				}
			default: throw new IllegalStateException();
			}
		}

	
	
	

	

	
	@Override
	protected int beforeVcf() {
		try {
			for(Path p: IOUtils.unrollPaths(this.databaseVCFPaths)) {
				this.databaseVCFS.add(new DatabaseVcf(p));
				}
			if(this.databaseVCFPaths.isEmpty()) {
				LOG.error("no database VCF was specified.");
				return -1;
				}
			
		}
		catch(final Throwable err ) {
			return -1;
			}
		return super.beforeVcf();
		}
	
	@Override
	protected void afterVcf() {
		
		for(DatabaseVcf v:this.databaseVCFS) {
			v.close();
			}
		super.afterVcf();
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator in, VariantContextWriter w) {
		try {
			final VCFHeader h0 = in.getHeader();
			final VCFHeader h = new VCFHeader(h0); 
			final VCFFilterHeaderLine filterHeaderLine;
			if(StringUtils.isBlank(this.filterName)) {
				filterHeaderLine = null;
				}
			else
				{
				filterHeaderLine = new VCFFilterHeaderLine(this.filterName, "Filtered "+getProgramName());
				h.addMetaDataLine(filterHeaderLine);
				}
			
			final VCFInfoHeaderLine countHdr = new VCFInfoHeaderLine("VCFIN_COUNT",1,VCFHeaderLineType.Integer,"The variant was found in 'x' other VCF using "+getProgramName());
			h.addMetaDataLine(countHdr);
			JVarkitVersion.getInstance().addMetaData(this, h);
			
			w.writeHeader(h);
			while(in.hasNext()) {
				final VariantContext ctx= in.next();
				final boolean is_sv1 = !StringUtils.isBlank(ctx.getAttributeAsString(VCFConstants.SVTYPE,""));
				
				int n_found=0;
				for(DatabaseVcf db:this.databaseVCFS) {
					boolean found = false;
					final String ctg = db.contigNameConverter.apply(ctx.getContig());
					if(!StringUtils.isBlank(ctg)) {
						try(CloseableIterator<VariantContext> iter = db.vcfReader.query(new Interval(ctg,ctx.getStart(),ctx.getEnd()))) {
							while(!found && iter.hasNext()) {
								final VariantContext vc2 = iter.next();
								if(!CoordMath.overlaps(ctx.getStart(), ctx.getEnd(), vc2.getStart(), vc2.getEnd())) continue;
								final boolean is_sv2 = !StringUtils.isBlank(vc2.getAttributeAsString(VCFConstants.SVTYPE,""));
								if(is_sv1 && is_sv2) {
									if(overlap_comparator.test(ctx.getStart(), ctx.getEnd(), vc2.getStart(), vc2.getEnd())) {
										found = true;
										}
									}
								else if(!is_sv1 && !is_sv2) {
									if(compare(this.matcher,ctx,vc2)) {
										
										}
									}
								else
									{
									/* both must have the same SV type */
									continue;
									}
								}
							
							}
						}
					if(found) n_found++;
 					}
				
				boolean keep = n_found>0;
				if(this.inverse) {
					keep= !keep;
					}
				
				// just skip the variant
				if(!keep && filterHeaderLine==null) continue;
				
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.attribute(countHdr.getID(), n_found);
				if(keep) {
					if(!ctx.isFiltered()) vcb.passFilters();
					}
				else
					{
					vcb.filter(filterHeaderLine.getID());
					}
				w.add(vcb.make());
				}
			return 0;
			}
		catch(final Throwable err ) {
			LOG.error(err);
			return -1;
		}
	}
		

	public static void main(final String[] args) {
		new VcfIn().instanceMainWithExit(args);
		}
	}
