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
package com.github.lindenb.jvarkit.net;

import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.igv.IgvConstants;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
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

private boolean isGrch37() {
	return hasDict() && SequenceDictionaryUtils.isGRCh37(this.dict);
	}
private boolean isGrch38() {
return hasDict() && SequenceDictionaryUtils.isGRCh38(this.dict);
}


public Set<LabelledUrl> of(final String id) {
	if(StringUtils.isBlank(id)) return Collections.emptySet();
	final Set<LabelledUrl> L = new LinkedHashSet<>();
	 _string(id,L);
	return L;
	}

public Set<LabelledUrl> of(final String columnName,final String id) {
	if(StringUtils.isBlank(columnName) || StringUtils.isBlank(id)) return Collections.emptySet();
	final Set<LabelledUrl> urls = new LinkedHashSet<>();
	if( columnName.equalsIgnoreCase("genename") ||
		columnName.equalsIgnoreCase("symbol")) {
		urls.add(new LabelledUrlImpl("NCBI gene","https://www.ncbi.nlm.nih.gov/gene/?term="+StringUtils.escapeHttp(id)));
		}
	else if(columnName.equalsIgnoreCase("Consequence") ||
			columnName.equalsIgnoreCase("so")) {
		for(final String ac:CharSplitter.of('&').split(id)) {
			if(StringUtils.isBlank(ac)) continue;
			if(ac.startsWith("SO:")) {
				urls.add(new LabelledUrlImpl("SO "+ac,"http://www.sequenceontology.org/browser/current_release/term/"+StringUtils.escapeHttp(ac)));
				} 
			else
				{
				urls.add(new LabelledUrlImpl("SO "+ac,"http://www.sequenceontology.org/browser/obob.cgi?rm=term_list&obo_query="+StringUtils.escapeHttp(ac)+"&release=current_svn"));
				}
			}
		}
	return urls;
	}


private void _string(final String str,final Set<LabelledUrl> urls) {
	if(StringUtils.isBlank(str)) return;
	if(this.rsIdPattern.matcher(str).matches())
		{
		urls.add(new LabelledUrlImpl("dbsnp","https://www.ncbi.nlm.nih.gov/snp/"+str.substring(2)));
		urls.add(new LabelledUrlImpl("opensnp","https://opensnp.org/snps/"+str.toLowerCase()));
		urls.add(new LabelledUrlImpl("Gnomad rs# GRCh37","https://gnomad.broadinstitute.org/variant/"+str.toLowerCase()+"?dataset=gnomad_r2_1"));
		urls.add(new LabelledUrlImpl("Gnomad rs# GRCh38","https://gnomad.broadinstitute.org/variant/"+str.toLowerCase()+"?dataset=gnomad_r3"));
		if(hasDict() && SequenceDictionaryUtils.isHuman(this.dict)) {
			urls.add(new LabelledUrlImpl("clinvar","https://www.ncbi.nlm.nih.gov/clinvar?term="+str.toLowerCase()+"%5BVariant%20ID%5D"));
			}
		}
	else if(this.ensemblPattern.matcher(str).matches())
		{
		urls.add(new LabelledUrlImpl("Ensembl","http://www.ensembl.org/Multi/Search/Results?species=all;idx=;q="+str.toUpperCase()+";species=;site=ensembl"));
		if(str.startsWith("ENSG")) {
			urls.add(new LabelledUrlImpl("Genbass GRCh38","https://genebass.org/gene/"+str+"?burdenSet=pLoF&phewasOpts=1&resultLayout=full"));
			}
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
	final String ensemblCtg = toEnsembl.apply(ctx.getContig());
	if(isGrch37() && !StringUtils.isBlank(ensemblCtg) && AcidNucleics.isATGC(ctx.getReference())) {
		for(final Allele alt: ctx.getAlternateAlleles()) {
			if(!AcidNucleics.isATGC(alt)) continue;
			//gnomad
			urls.add(new LabelledUrlImpl("Variant Gnomad 2.1 " + alt.getDisplayString(),"https://gnomad.broadinstitute.org/variant/"+
				StringUtils.escapeHttp(ensemblCtg) + "-" + ctx.getStart() +"-"+ctx.getReference().getDisplayString()+"-"+alt.getDisplayString()+"?dataset=gnomad_r2_1"
				));
			
			
			// Bibliome
			urls.add(new LabelledUrlImpl("Bibliome " + alt.getDisplayString(),"https://bibliome.ai/variant/"+
				ensemblCtg +
				"-"+ctx.getStart()+
				"-"+
				ctx.getReference().getDisplayString()+
				"-"+
				alt.getDisplayString()
				));
			// marvel https://twitter.com/julawang/status/1094666160711323649
			urls.add(new LabelledUrlImpl("Marrvel "+ alt.getDisplayString(),"http://marrvel.org/search/variant/"+
				ensemblCtg +
				"-"+ctx.getStart()+
				StringUtils.escapeHttp("+")+
				ctx.getReference().getDisplayString()+
				StringUtils.escapeHttp(">")+
				alt.getDisplayString()
				));
			}
		
		
		
		}
	if(isGrch38() && ! StringUtils.isBlank(ensemblCtg) && AcidNucleics.isATGC(ctx.getReference())) {
		for(final Allele alt: ctx.getAlternateAlleles()) {
			if(!AcidNucleics.isATGC(alt)) continue;
			urls.add(new LabelledUrlImpl("Variant Gnomad 3 " + alt.getDisplayString(),"https://gnomad.broadinstitute.org/variant/"+
				StringUtils.escapeHttp(ensemblCtg) + "-" + ctx.getStart() +"-"+ctx.getReference().getDisplayString()+"-"+alt.getDisplayString()+"?dataset=gnomad_r3"
				));
			}
		}
	
	
	//beacon , //varsome
	for(int side=0;side<2 && !StringUtils.isBlank(ensemblCtg);++side)
		{
		if(side==0 && !isGrch37()) continue;
		if(side==1 && !isGrch38()) continue;
		for(final Allele alt: ctx.getAlternateAlleles())
			{
			if(ctx.getReference().isSymbolic() || alt.isSymbolic()) continue;
			//https://beacon-network.org/#/search?pos=114267128&chrom=4&allele=A&ref=G&rs=GRCh37
			urls.add(new LabelledUrlImpl("Beacon " + alt.getDisplayString(),
					"https://beacon-network.org/#/search?chrom="+
					StringUtils.escapeHttp(ensemblCtg) +
					"&pos="+ctx.getStart()+
					"&ref="+ ctx.getReference().getDisplayString()+
					"&allele="+ alt.getDisplayString()+
					"&rs="+ (side==0?"GRCh37":"GRCh38")
					));
			urls.add(new LabelledUrlImpl("Varsome " + alt.getDisplayString(),
					"https://varsome.com/variant/"+
					(side==0?"hg19/":"hg38/")+
					StringUtils.escapeHttp(ensemblCtg) + "-"+
					ctx.getStart()+ "-"+
					ctx.getReference().getDisplayString()+"-"+
					alt.getDisplayString()
					));

			}
		}
		
	_locatable(ctx, urls);
	}

private void _sam(final SAMRecord rec,final Set<LabelledUrl> urls) {
	if(rec==null) return;
	if(!rec.getReadUnmappedFlag()) {
		_locatable(rec, urls);
		}
	}

private void _interval(final Locatable loc,final Set<LabelledUrl> urls) {
	if(loc==null) return;
	int extend = 100;
	final int xstart1 = Math.max(loc.getStart()-extend,1);
	final int xend1 = loc.getEnd()+1;
	urls.add(new LabelledUrlImpl("IGV","https://"+ IgvConstants.DEFAULT_HOST +":"+IgvConstants.DEFAULT_PORT + "/goto?locus="+
			StringUtils.escapeHttp(loc.getContig()) + "%3A" +xstart1 +"-"+loc.getEnd()
			));
	
	final String ensemblCtg = toEnsembl.apply(loc.getContig());
	if(isGrch37() && ! StringUtils.isBlank(ensemblCtg)) {
		urls.add(new LabelledUrlImpl("Region Gnomad 2.1","https://gnomad.broadinstitute.org/region/"+
			StringUtils.escapeHttp(ensemblCtg) + "-" + xstart1 +"-"+ xend1 +"?dataset=gnomad_r2_1"
			));
		
		urls.add(new LabelledUrlImpl("clinvar 37","https://www.ncbi.nlm.nih.gov/clinvar/?term="+
				ensemblCtg +
				"%5Bchr%5D+AND+"+ loc.getStart()+"%3A"+ loc.getEnd()+"%5Bchrpos37%5D"
				));
			}
	if(isGrch38() && ! StringUtils.isBlank(ensemblCtg)) {
		urls.add(new LabelledUrlImpl("Region Gnomad 3","https://gnomad.broadinstitute.org/region/"+
			StringUtils.escapeHttp(ensemblCtg) + "-" + xstart1 +"-"+ xend1 +"?dataset=gnomad_3"
			));
		}
	
	final String ucscCtg =  toUcsc.apply(loc.getContig());
	if(isGrch37() && ! StringUtils.isBlank(ucscCtg)) {
		urls.add(new LabelledUrlImpl("UCSC hg19","http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg19&position="+
			StringUtils.escapeHttp(ucscCtg) + "%3A" + xstart1 +"-"+ xend1 + "&highlight="+
			StringUtils.escapeHttp(ucscCtg) + "%3A" + loc.getStart() +"-"+ loc.getEnd()
			));
		
		if(loc.getLengthOnReference()>1) {
			
			urls.add(new LabelledUrlImpl("dgv","http://dgv.tcag.ca/gb2/gbrowse/dgv2_hg19?name="+
					StringUtils.escapeHttp(ucscCtg) + 
					"%3A"+loc.getStart()+"-"+loc.getEnd() + ";search=Search"
					));
			}
		
		}
	if(isGrch38() && ! StringUtils.isBlank(ucscCtg)) {
		urls.add(new LabelledUrlImpl("UCSC hg38","http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg38&position="+
			StringUtils.escapeHttp(ucscCtg) + "%3A" + xstart1 +"-"+ xend1 + "&highlight="+
					StringUtils.escapeHttp(ucscCtg) + "%3A" + loc.getStart() +"-"+ loc.getEnd()
			));
		if(loc.getLengthOnReference()>1) {
			urls.add(new LabelledUrlImpl("decipher","https://decipher.sanger.ac.uk/search?q="+
			StringUtils.escapeHttp(ucscCtg) + 
			"%3A"+loc.getStart()+"-"+loc.getEnd()));
			}
		}
	
	if(loc.getLengthOnReference()>1) {
		for(int i=0;i< 2;++i) {
			if(i==0 && !isGrch37()) continue;
			if(i==1 && !isGrch38()) continue;
			urls.add(new LabelledUrlImpl("Hi-C","http://promoter.bx.psu.edu/hi-c/view.php?method=Hi-C&species=human&assembly="+(i==0?"hg19":"hg38")+
					"&source=inside&tissue=GM12878&type=Lieberman-raw&resolution=25&c_url=&transfer=&gene=&chr="+
					StringUtils.escapeHttp(loc.getContig())+"&start="+loc.getStart()+"&end="+loc.getEnd()+"&sessionID=&browser=none"
					));		
			}
		}
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
