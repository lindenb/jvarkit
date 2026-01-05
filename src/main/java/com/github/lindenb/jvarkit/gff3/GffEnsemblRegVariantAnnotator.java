/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.gff3;

import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.github.lindenb.jvarkit.iterator.LineIterators;
import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.tabix.AbstractTabixVariantAnnotator;

import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/*** tabix annotator using Ensembl GFF to annotate SV (e/g: https://ftp.ensembl.org/pub/grch37/current/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff.gz ) */
public class GffEnsemblRegVariantAnnotator extends AbstractTabixVariantAnnotator {
	private static final Logger LOG = Logger.of(GffEnsemblRegVariantAnnotator.class);
	public static final String OPT_DESC="GFF file from Ensembl, indexed with Tabix. e:g https://ftp.ensembl.org/pub/grch37/current/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff.gz";
	private final Set<String> undefined_types = new HashSet<>();
	private final Map<String,VCFInfoHeaderLine> type2Flag = new HashMap<>();
	private final Gff3Codec gff3Codec;
	public GffEnsemblRegVariantAnnotator(final String uri, AttributeMap meta) throws IOException {
		super(uri);
		if(meta==null) meta = AttributeMap.empty();
		this.gff3Codec = new Gff3Codec(Gff3Codec.DecodeDepth.SHALLOW);
		}

	@Override
	public void fillHeader(final VCFHeader header) {
		if(super.tabixReader==null) return;
		String helpCmd = " in " + getUri();
		final String[] types = new String[] {
				"CTCF_binding_site", "enhancer", "open_chromatin_region",
				"promoter", "promoter_flanking_region", "TF_binding_site"
				};
		for(String type: types) {
			final VCFInfoHeaderLine h = new VCFInfoHeaderLine(
					"IN_" + type.toUpperCase(),1,
					VCFHeaderLineType.Flag, "variant was found overlaping a  "+ type+ helpCmd
					);
			this.type2Flag.put(type,h);
			header.addMetaDataLine(h);
			}
		}
	@Override
	public List<VariantContext> annotate(final VariantContext ctx) throws IOException  {
		if(super.tabixReader ==null || !hasContig(ctx)) return Collections.singletonList(ctx);
		final Set<String> all_types = new HashSet<>();
		final TabixReader.Iterator r = this.tabixReader.query(contig(ctx), ctx.getStart(), ctx.getEnd());
		for(;;) {
			final String line = r.next();
			if(line==null) break;
			final Gff3Feature feat = this.gff3Codec.decode(LineIterators.of(Collections.singletonList(line)));
			if(feat==null) continue;
			all_types.add(feat.getType());
			}
			
		if(all_types.isEmpty()) {
			return Collections.singletonList(ctx);
			}
		
		final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
		for(final String type : all_types) {
			if(!this.type2Flag.containsKey(type)) {
				if(undefined_types.add(type)) {
					LOG.warning("undefined reg type "+type+" in "+getUri());
					}
				continue;
				}
			vcb.attribute(this.type2Flag.get(type).getID(), true);
			}
		return Collections.singletonList(vcb.make());
		}
	}
