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
package com.github.lindenb.jvarkit.tools.iranome;


import java.io.BufferedReader;
import java.io.InputStream;
import java.util.Arrays;
import java.util.Map;
import java.util.Optional;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.vcf.BufferedVCFReader;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFIterator;
import mjson.Json;

public class IranomeScrapper extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.build(IranomeScrapper.class).make();

	private BufferedVCFReader database;
	private VariantContextWriter saveDatabaseWriter = null;
	private static final String NO_FOUND_IN_IRANOME_FILTER = "NOT_FOUND_IN_IRANOME";
	
	private class VarInfo {
		int ac;
		int an;
		Double af;
		}
	
	private String toIranomeContig(final String c) {
		if(c.startsWith("chr")) return c.substring(3);
		return c;
		}
	
	private Optional<VarInfo> getIranomeVariant(final VariantContext userCtx,final Allele alt) {
		final String iranomeCtg = toIranomeContig(userCtx.getContig());
		if(StringUtils.isBlank(iranomeCtg)) return null;
		Locatable userLoc = new SimpleInterval(iranomeCtg, userCtx.getStart(), userCtx.getEnd());
		try(CloseableIterator<VariantContext> iter=this.database.query(userLoc)) {
			while(iter.hasNext()) {
				final VariantContext ctx = iter.next();
				if(!ctx.getReference().equals(userCtx.getReference())) continue;
				if(!ctx.hasAlternateAllele(alt)) continue;
				if(ctx.getFilters().contains(NO_FOUND_IN_IRANOME_FILTER)) {
					return Optional.empty();
					}
				}
			}
		final VariantContextBuilder vcb = new VariantContextBuilder(null, iranomeCtg, userCtx.getStart(), userCtx.getEnd(), Arrays.asList(userCtx.getReference(),alt));
		vcb.attribute("date", StringUtils.now());

		int got_data = 0;
		String url = "http://www.iranome.ir/variant/"+iranomeCtg+"-"+userCtx.getStart()+"-"+userCtx.getReference().getDisplayString()+"-"+alt.getDisplayString();
		try(InputStream is = openURL(url)) {
			final String window_variant= "window.variant = ";
			final String end_json = "};";
			final String html = IOUtil.slurp(is);
			int i1= html.indexOf(window_variant);
			final VarInfo varInfo = new VarInfo();
			if(i1!=-1) {
				int i2 = html.indexOf(end_json);
				if(i1==-1) return null;
				final String jsonString = html.substring(i1+window_variant.length(),i2+1);
				mjson.Json parsed = mjson.Json.read(jsonString);
				final Map<String,Json> as_map = parsed.asJsonMap();
				if(as_map.containsKey("allele_count")) {
					varInfo.ac = as_map.get("allele_count").asInteger();
					vcb.attribute(VCFConstants.ALLELE_COUNT_KEY,varInfo.ac);
					}
				if(as_map.containsKey("allele_count")) {
					varInfo.an = as_map.get("allele_count").asInteger();
					vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY,varInfo.an);
					}
				
				if(as_map.containsKey("allele_freq")) {
					varInfo.af = as_map.get("allele_freq").asDouble();
					vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY,varInfo.af);
					}
				}
			saveDatabaseWriter.add(vcb.make());
			return Optional.of(varInfo);
			}
		vcb.filter(NO_FOUND_IN_IRANOME_FILTER);
		saveDatabaseWriter.add(vcb.make());
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin, VariantContextWriter out) {
		try {
			
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	public static void main(String[] args) {
		new IranomeScrapper().instanceMainWithExit(args);
		}
	}
