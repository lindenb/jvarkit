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
package com.github.lindenb.jvarkit.net;

import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.igv.IgvConstants;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * URL supplier
 */
public class UrlSupplier {
private final SAMSequenceDictionary dict;
private final ContigNameConverter toUcsc = ContigNameConverter.createConvertToUcsc();
private final ContigNameConverter toEnsembl = ContigNameConverter.createConvertToEnsembl();
private final Pattern rsIdPattern = Pattern.compile("[rR][sS][0-9]+");
private final Pattern ensemblPattern = Pattern.compile("ENS(EST[TGP]|[TPGR])[0-9]+");
private final Pattern ccdsPattern = Pattern.compile("CCDS[0-9\\.]+");
private final Pattern ncbiNucPattern = Pattern.compile("[XN][MR]_[0-9\\.]+");
private final Pattern ncbiProtPattern = Pattern.compile("[NX]P_[0-9\\.]+");


public static interface LabelledUrl
	{
	public String getLabel();
	public String getUrl();
	}

public UrlSupplier(final SAMSequenceDictionary dict) {
	this.dict = dict;
	}

public UrlSupplier() {
	this(new SAMSequenceDictionary());
	}
private boolean hasDict() {
	return this.dict!=null && !this.dict.isEmpty();
}

public Set<LabelledUrl> of(final String id) {
	if(StringUtils.isBlank(id)) return Collections.emptySet();
	final Set<LabelledUrl> L = new LinkedHashSet<>();
	 _string(id,L);
	return L;
	}

private void _string(final String str,final Set<LabelledUrl> urls) {
	if(StringUtils.isBlank(str)) return;
	if(this.rsIdPattern.matcher(str).matches())
		{
		urls.add(new LabelledUrlImpl("dbsnp","https://www.ncbi.nlm.nih.gov/snp/"+str.substring(2)));
		urls.add(new LabelledUrlImpl("opensnp","https://opensnp.org/snps/"+str.toLowerCase()));
		if(hasDict() && SequenceDictionaryUtils.isHuman(this.dict)) {
			urls.add(new LabelledUrlImpl("clinvar","https://www.ncbi.nlm.nih.gov/clinvar?term="+str.toLowerCase()+"%5BVariant%20ID%5D"));
			}
		}
	else if(this.ensemblPattern.matcher(str).matches())
		{
		urls.add(new LabelledUrlImpl("Ensembl","http://www.ensembl.org/Multi/Search/Results?species=all;idx=;q="+str.toUpperCase()+";species=;site=ensembl"));
		}
	else if(this.ccdsPattern.matcher(str).matches()) {
		urls.add(new LabelledUrlImpl("CCDS","https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&GO=MainBrowse&DATA="+str));
		}
	else if(this.ncbiNucPattern.matcher(str).matches()) {
		urls.add(new LabelledUrlImpl("NCBI","https://www.ncbi.nlm.nih.gov/nuccore/"+str));
		}
	else if(this.ncbiProtPattern.matcher(str).matches()) {
		urls.add(new LabelledUrlImpl("NCBI","https://www.ncbi.nlm.nih.gov/protein/"+str));
		}
	else if(IOUtil.isUrl(str))
		{
		urls.add(new LabelledUrlImpl("url",str));
		}
	}


private void _variant(final VariantContext ctx,final Set<LabelledUrl> urls) {
	if(ctx==null) return;
	if(ctx.hasID()) {
		_string(ctx.getID(),urls);
		}
	_locatable(ctx, urls);
	}

private void _sam(final SAMRecord rec,final Set<LabelledUrl> urls) {
	if(rec==null) return;
	_locatable(rec, urls);
	}

private void _interval(final Locatable loc,final Set<LabelledUrl> urls) {
	if(loc==null) return;
	urls.add(new LabelledUrlImpl("IGV","https://"+ IgvConstants.DEFAULT_HOST +":"+IgvConstants.DEFAULT_PORT + "/goto?locus="+
			StringUtils.escapeHttp(loc.getContig()) + "%3A" + loc.getStart() +"-"+loc.getEnd()
			));
	}


private void _locatable(final Locatable loc,final Set<LabelledUrl> urls) {
	if(loc==null) return;
	_interval(loc,urls);
	}

public Set<LabelledUrl> of(final Locatable loc) {
	if(loc==null) return Collections.emptySet();
	final Set<LabelledUrl> urls = new LinkedHashSet<>();
	
	if(loc instanceof VariantContext) {
		_variant(VariantContext.class.cast(loc), urls);
		}
	else if(loc instanceof SAMRecord) {
		_sam(SAMRecord.class.cast(loc), urls);
		}
	else
		{
		_locatable(loc,urls);
		}
	return urls;
	}

private static class LabelledUrlImpl implements LabelledUrl
	{
	private final String label;
	private final String url;
	LabelledUrlImpl(final String label,final String url) {
		this.label = label;
		this.url = url;
		}
	public String getLabel() { return label;}
	public String getUrl() { return url;}
	@Override
	public int hashCode() {
		return getUrl().hashCode();
		}
	@Override
	public boolean equals(final Object obj) {
		if(obj==null || !(obj instanceof LabelledUrl)) return false;
		if(obj==this) return true;
		return getUrl().equals(LabelledUrl.class.cast(obj).getUrl());
		}
	@Override
	public String toString() {
		return getUrl();
		}
	}
}
