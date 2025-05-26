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
package com.github.lindenb.jvarkit.net;

import java.net.URL;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.igv.IgvConstants;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.LocatableUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;

/**
 * URL supplier
 */
public class UrlSupplier {
private static final String GNOMAD_HG38 = "gnomad_r4";

private final SAMSequenceDictionary dict;
private final ContigNameConverter toUcsc = ContigNameConverter.createConvertToUcsc();
private final ContigNameConverter toEnsembl = ContigNameConverter.createConvertToEnsembl();
private final Pattern rsIdPattern = Pattern.compile("[rR][sS][0-9]+");
private final Pattern ensemblPattern = Pattern.compile("ENS[TGP][0-9]+(\\.[0-9]+)?");
private final Pattern ccdsPattern = Pattern.compile("CCDS[0-9\\.]+");
private final Pattern ncbiNucPattern = Pattern.compile("[XN][MR]_[0-9\\.]+");
private final Pattern ncbiProtPattern = Pattern.compile("[NX]P_[0-9\\.]+");
private final Pattern hgncPattern = Pattern.compile("[hH][gG][nN][cC]:[0-9]+");
public static interface LabelledUrl extends Comparable<LabelledUrl>
	{
	/** short name for this url : database... */
	public String getLabel();
	/** target for this url: gene Name, gene id, interval... */
	public String getTarget();
	/** the url itself */
	public String getUrl();
	/** return domain of getURL() */
	public default String getDomain() {
		final String s= getUrl();
		try {
			final java.net.URL u = new URL(s);
			return u.getProtocol()+"://"+u.getHost();
			}
		catch(final Throwable err) {
			return s;
			}
		}
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
		columnName.equalsIgnoreCase("gene_name") ||
		columnName.equalsIgnoreCase("symbol")) {
		urls.add(new LabelledUrlImpl("Pharos",id,"https://pharos.nih.gov/diseases?associatedTarget="+StringUtils.escapeHttp(id)));
		urls.add(new LabelledUrlImpl("NCBI gene",id,"https://www.ncbi.nlm.nih.gov/gene/?term="+StringUtils.escapeHttp(id)));
		urls.add(new LabelledUrlImpl("OMIM",id,"https://www.omim.org/search?search="+StringUtils.escapeHttp(id)));
		urls.add(new LabelledUrlImpl("hugeamp",id,"https://cvd.hugeamp.org/gene.html?gene="+StringUtils.escapeHttp(id)));
		urls.add(new LabelledUrlImpl("Archs4",id,"https://maayanlab.cloud/archs4/gene/"+StringUtils.escapeHttp(id)));
		urls.add(new LabelledUrlImpl("Enrichr",id,"https://maayanlab.cloud/Enrichr/#find!gene="+StringUtils.escapeHttp(id)));
		urls.add(new LabelledUrlImpl("Biogps",id,"http://biogps.org/#goto=search=&query="+StringUtils.escapeHttp(id)));
		urls.add(new LabelledUrlImpl("Gene ResearchAllOfUs",id,"https://databrowser.researchallofus.org/genomic-variants/"+ StringUtils.escapeHttp(id)));
		// intron retention associated variants
		urls.add(new LabelledUrlImpl("IRAVs",id,"https://iravdb.io/gene/"+ StringUtils.escapeHttp(id)));
		//finngen
		urls.add(new LabelledUrlImpl("Finngen",id,"https://public-metaresults-fg-ukbb.finngen.fi/gene/"+ StringUtils.escapeHttp(id)));
		
		urls.add(new LabelledUrlImpl("hugeamp",id,"https://hugeamp.org:8000/research.html?ancestry=mixed&cohort=AoU_250k&file=600Traits.csv&gene="+ StringUtils.escapeHttp(id) +"&pageid=600_traits_app&phenotype=phecode_425.0"));
		urls.add(new LabelledUrlImpl("intogen",id,"https://www.intogen.org/search?gene="+ StringUtils.escapeHttp(id)));
		urls.add(new LabelledUrlImpl("ComplexPortal",id,"https://www.ebi.ac.uk/complexportal/complex/search?query="+ StringUtils.escapeHttp(id)));

		 if(isGrch38()) {
			urls.add(new LabelledUrlImpl("TogoVar",id,"https://grch38.togovar.org/?mode=simple&term="+ StringUtils.escapeHttp(id) ));
			}
		 // Varaico: Literature-Extracted Genetic Variants
		 if(isGrch38() || isGrch37()) {
			urls.add(new LabelledUrlImpl("Varaico",id,
					"https://varaico.com/?gene="+ StringUtils.escapeHttp(id) +"&assembly="+ (isGrch38()?"hg38":(isGrch37()?"hg19":"undefined"))));
			}
		 
		 
		 
		}
	else if( columnName.equalsIgnoreCase("hgnc") && (StringUtils.isInteger(id)  || id.toUpperCase().startsWith("HGNC:"))) {
		final String hgnc = (StringUtils.isInteger(id)?"HGNC:":"") + id.toUpperCase();
		urls.add(new LabelledUrlImpl("GenCC",id,"https://search.thegencc.org/genes/"+hgnc));
		urls.add(new LabelledUrlImpl("Monarch",id,"https://monarchinitiative.org/gene/" + hgnc));
		}
	else if( columnName.equalsIgnoreCase("omim")) {
		urls.add(new LabelledUrlImpl("OMIM",id,"https://www.omim.org/entry/"+id));
		}
	else if(columnName.equalsIgnoreCase("Consequence") ||
			columnName.equalsIgnoreCase("so")) {
		for(final String ac:CharSplitter.of('&').split(id)) {
			if(StringUtils.isBlank(ac)) continue;
			if(ac.startsWith("SO:")) {
				urls.add(new LabelledUrlImpl("SO "+ac,id,"http://www.sequenceontology.org/browser/current_release/term/"+StringUtils.escapeHttp(ac)));
				} 
			else
				{
				urls.add(new LabelledUrlImpl("SO "+ac,id,"http://www.sequenceontology.org/browser/obob.cgi?rm=term_list&obo_query="+StringUtils.escapeHttp(ac)+"&release=current_svn"));
				}
			}
		}
	else if(columnName.equalsIgnoreCase("UniProtKB")) {
		urls.add(new LabelledUrlImpl("functionome",id,"https://functionome.geneontology.org/gene/UniProtKB:"+id));
		urls.add(new LabelledUrlImpl("amigo",id,"https://amigo.geneontology.org/amigo/gene_product/UniProtKB:"+id));
		}
	return urls;
	}

private void _gff3(final Gff3Feature feat,final Set<LabelledUrl> urls) {
	feat.getAttribute("gene_id").stream().forEach(X->_string(X,urls));;
	feat.getAttribute("transcript_id").stream().forEach(X->_string(X,urls));;
	feat.getAttribute("protein_id").stream().forEach(X->_string(X,urls));;
	feat.getAttribute("gene_name").stream().forEach(X->urls.addAll(of("gene_name",X)));;
	feat.getAttribute("CCDS").stream().forEach(X->urls.addAll(of("CCDS",X)));;
	_locatable(feat, urls);
	}
private void _string(final String str,final Set<LabelledUrl> urls) {
	if(StringUtils.isBlank(str)) return;
	else if(str.startsWith("UniProtKB:")) {
		urls.addAll(of("UniProtKB",str.substring(10)));
		}
	else if(this.hgncPattern.matcher(str).matches()) {
		urls.addAll(of("hgnc",str));
		}
	if(this.rsIdPattern.matcher(str).matches())
		{
		urls.add(new LabelledUrlImpl("dbsnp",str,"https://www.ncbi.nlm.nih.gov/snp/"+str.substring(2)));
		urls.add(new LabelledUrlImpl("Japan-Omics",str,"https://japan-omics.jp/variant?input_value="+str.substring(2)));
		urls.add(new LabelledUrlImpl("opensnp",str,"https://opensnp.org/snps/"+str.toLowerCase()));
		urls.add(new LabelledUrlImpl("Gnomad rs# GRCh37",str,"https://gnomad.broadinstitute.org/variant/"+str.toLowerCase()+"?dataset=gnomad_r2_1"));
		urls.add(new LabelledUrlImpl("Gnomad rs# GRCh38",str,"https://gnomad.broadinstitute.org/variant/"+str.toLowerCase()+"?dataset="+GNOMAD_HG38));
		if(hasDict() && SequenceDictionaryUtils.isHuman(this.dict)) {
			urls.add(new LabelledUrlImpl("clinvar",str,"https://www.ncbi.nlm.nih.gov/clinvar?term="+str.toLowerCase()+"%5BVariant%20ID%5D"));
			}
		if(isGrch38()) {
			urls.add(new LabelledUrlImpl("TogoVar",str,"https://grch38.togovar.org/?mode=simple&term="+ str ));
			}
		}
	else if(this.ensemblPattern.matcher(str).matches())
		{
		urls.add(new LabelledUrlImpl("Ensembl",str,"http://www.ensembl.org/Multi/Search/Results?species=all;idx=;q="+str.toUpperCase()+";species=;site=ensembl"));
		if(str.startsWith("ENSG")) {
			urls.add(new LabelledUrlImpl("OpenTargets",str,"https://genetics.opentargets.org/gene/"+str));
			urls.add(new LabelledUrlImpl("Genbass GRCh38",str,"https://genebass.org/gene/"+str+"?burdenSet=pLoF&phewasOpts=1&resultLayout=full"));
			urls.add(new LabelledUrlImpl("Protein Atlas",str,"https://www.proteinatlas.org/"+str));
			urls.add(new LabelledUrlImpl("Japan-Omics",str,"https://japan-omics.jp/gene/JCTF?input_value="+str));
			}
		}
	else if(this.ccdsPattern.matcher(str).matches()) {
		urls.add(new LabelledUrlImpl("CCDS",str,"https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&GO=MainBrowse&DATA="+str));
		}
	else if(this.ncbiNucPattern.matcher(str).matches()) {
		urls.add(new LabelledUrlImpl("NCBI",str,"https://www.ncbi.nlm.nih.gov/nuccore/"+str));
		}
	else if(this.ncbiProtPattern.matcher(str).matches()) {
		urls.add(new LabelledUrlImpl("NCBI",str,"https://www.ncbi.nlm.nih.gov/protein/"+str));
		}
	else if(str.matches("nsv[0-9]+")) {
		if(isGrch37()) {
			urls.add(new LabelledUrlImpl("DGV",str,"http://dgv.tcag.ca/gb2/gbrowse/dgv2_hg19?name="+str+"&search=Search"));
			}
		else if(isGrch38()) {
			urls.add(new LabelledUrlImpl("DGV",str,"http://dgv.tcag.ca/gb2/gbrowse/dgv2_hg38?name="+str+"&search=Search"));
			}
		}
	else if(str.matches("tgv[0-9]+")) {
		if(isGrch38()) {
			urls.add(new LabelledUrlImpl("TogoVar",str,"https://grch38.togovar.org/variant/"+ str ));
			}
		}
	else if(IOUtil.isUrl(str))
		{
		urls.add(new LabelledUrlImpl("url",str,str));
		}
	}

/* TODO : https://japan-omics.jp/ */

private void _variant(final VariantContext ctx,final Set<LabelledUrl> urls) {
	if(ctx==null) return;
	if(ctx.hasID()) {
		for(final String id: CharSplitter.COMMA.split(ctx.getID())) {
			_string(id,urls);
			}
		}
	final String ensemblCtg = toEnsembl.apply(ctx.getContig());
	final String ucscCtg = toUcsc.apply(ctx.getContig());
	final String variantid = ctx.getContig()+":"+ctx.getStart()+":"+ctx.getReference().getDisplayString();
	
	// FREX
	if(isGrch37() && ctx.getStart()==ctx.getEnd()) {
		urls.add(new LabelledUrlImpl("FreX",
				variantid,
				"http://lysine.univ-brest.fr/FrExAC/showAnnotations?chr="+
				StringUtils.escapeHttp(toUcsc.apply(ctx.getContig())) +
				"&position=" + ctx.getStart()
				));
		}
	
	// Franklin
	if(isGrch37()) {
		if(ctx.hasAttribute(VCFConstants.SVTYPE)) {
				final String svType = ctx.getAttributeAsString(VCFConstants.SVTYPE, "SV");
				urls.add(new LabelledUrlImpl("genoox",
						variantid,
						"https://franklin.genoox.com/clinical-db/variant/sv/"+
						StringUtils.escapeHttp(ctx.getContig()) + "-" + ctx.getStart()+"-"+ctx.getEnd()+"-"+svType
						));
			} else if(AcidNucleics.isATGC(ctx.getReference())) {
				for(final Allele alt: ctx.getAlternateAlleles()) {
					if(!AcidNucleics.isATGC(alt)) continue;
					urls.add(new LabelledUrlImpl("genoox",
							variantid,
							"https://franklin.genoox.com/clinical-db/variant/snp/"+
							StringUtils.escapeHttp(ctx.getContig()) + "-" + ctx.getStart()+"-"+ctx.getEnd()+"-"+ctx.getReference().getDisplayString()+"-"+alt.getDisplayString()
							));
				}
			}
		}
	
	// popgen
	if(isGrch37()) {
		urls.add(new LabelledUrlImpl("popgen.uchicago.edu",
				variantid,
				"https://popgen.uchicago.edu/ggv/?data=%221000genomes%22&chr="+
				StringUtils.escapeHttp(ctx.getContig()) + "&pos=" + ctx.getStart()
				));
		}
	
	if(isGrch37() && !StringUtils.isBlank(ensemblCtg) && AcidNucleics.isATGC(ctx.getReference())) {
		for(final Allele alt: ctx.getAlternateAlleles()) {
			if(!AcidNucleics.isATGC(alt)) continue;
			//gnomad
			urls.add(new LabelledUrlImpl("Variant Gnomad 2.1 " + alt.getDisplayString(),variantid+"/"+alt.getDisplayString(),"https://gnomad.broadinstitute.org/variant/"+
				StringUtils.escapeHttp(ensemblCtg) + "-" + ctx.getStart() +"-"+ctx.getReference().getDisplayString()+"-"+alt.getDisplayString()+"?dataset=gnomad_r2_1"
				));
			
			
			// Bibliome
			urls.add(new LabelledUrlImpl("Bibliome " + alt.getDisplayString(),
				variantid+"/"+alt.getDisplayString(),
				"https://bibliome.ai/variant/"+
				ensemblCtg +
				"-"+ctx.getStart()+
				"-"+
				ctx.getReference().getDisplayString()+
				"-"+
				alt.getDisplayString()
				));
			// marvel https://twitter.com/julawang/status/1094666160711323649
			urls.add(new LabelledUrlImpl("Marrvel "+ alt.getDisplayString(),
				variantid+"/"+alt.getDisplayString(),
				"http://marrvel.org/search/variant/"+
				ensemblCtg +
				"-"+ctx.getStart()+
				StringUtils.escapeHttp("+")+
				ctx.getReference().getDisplayString()+
				StringUtils.escapeHttp(">")+
				alt.getDisplayString()
				));
			}
		}
	if(isGrch38() && AcidNucleics.isATGC(ctx.getReference())) {
		for(final Allele alt: ctx.getAlternateAlleles()) {
			if(!AcidNucleics.isATGC(alt)) continue;
			if(!StringUtils.isBlank(ensemblCtg)) {
				urls.add(new LabelledUrlImpl("Variant Gnomad 3 " + alt.getDisplayString(),
						variantid+"/"+alt.getDisplayString(),
						"https://gnomad.broadinstitute.org/variant/"+
						StringUtils.escapeHttp(ensemblCtg) + "-" + ctx.getStart() +"-"+ctx.getReference().getDisplayString()+"-"+alt.getDisplayString()+"?dataset=" +GNOMAD_HG38
						));
				urls.add(new LabelledUrlImpl("OpenTargets " + alt.getDisplayString(),
						variantid+"/"+alt.getDisplayString(),
						"https://genetics.opentargets.org/variant/"+
						String.join("_",StringUtils.escapeHttp(ensemblCtg) ,""+ctx.getStart(),ctx.getReference().getDisplayString(),alt.getDisplayString())
						));
				
				urls.add(new LabelledUrlImpl("MexicoCPS" + alt.getDisplayString(),
						variantid+"/"+alt.getDisplayString(),
						"https://rgc-mcps.regeneron.com/variant/"+
						StringUtils.escapeHttp(ensemblCtg) + ":" + ctx.getStart() +":"+ctx.getReference().getDisplayString()+":"+alt.getDisplayString()
						));
				
				urls.add(new LabelledUrlImpl("rgc-research" + alt.getDisplayString(),
						variantid+"/"+alt.getDisplayString(),
						"https://rgc-research.regeneron.com/me/variant/"+
						StringUtils.escapeHttp(ensemblCtg) + ":" + ctx.getStart() +":"+ctx.getReference().getDisplayString()+":"+alt.getDisplayString()
						));
				}
			
			if(!StringUtils.isBlank(ucscCtg)) {
				urls.add(new LabelledUrlImpl("Decaf Decode",
						variantid+"/"+alt.getDisplayString(),
						"https://decaf.decode.com/variant/"+
						StringUtils.escapeHttp(ucscCtg) +":"+ctx.getStart()+":SG"
						));

				
				urls.add(new LabelledUrlImpl("Japan-Omics",
						variantid+"/"+alt.getDisplayString(),
						"https://japan-omics.jp/variant?input_value="+
						StringUtils.escapeHttp(
								ucscCtg+":"+
								ctx.getStart()+":"+
								ctx.getReference().getBaseString()+":"+
								alt.getBaseString()
								) 
						));
				
				}
			urls.add(new LabelledUrlImpl("Variant ResearchAllOfUs",
					variantid+"/"+alt.getDisplayString(),
					"https://databrowser.researchallofus.org/genomic-variants/"+
					StringUtils.escapeHttp(ctx.getContig()) + "-" + ctx.getStart() +"-"+ctx.getReference().getDisplayString()+"-"+alt.getDisplayString()
					));


			if(!StringUtils.isBlank(ucscCtg)) {
				urls.add(new LabelledUrlImpl("AF.ukbiobank",
						variantid+"/"+alt.getDisplayString(),
						"https://afb.ukbiobank.ac.uk/variant/"+
						StringUtils.escapeHttp(ucscCtg) + "-" + ctx.getStart() +"-"+ctx.getReference().getDisplayString()+"-"+alt.getDisplayString()
						));
				}
			
			urls.add(new LabelledUrlImpl("Genebe",
					variantid+"/"+alt.getDisplayString(),
					"https://genebe.net/variant/hg38/"+
					StringUtils.escapeHttp(ctx.getContig()) + ":" +
					ctx.getStart() +"-"+ctx.getReference().getDisplayString()+">"+alt.getDisplayString()
					));
			}
		}
	
	//spliceAI
	if(AcidNucleics.isATGC(ctx.getReference()) && (isGrch37() || isGrch38())) {
		for(final Allele alt: ctx.getAlternateAlleles()) {
			if(!AcidNucleics.isATGC(alt)) continue;
			urls.add(new LabelledUrlImpl("SpliceAI",
					variantid+"/"+alt.getDisplayString(),
					"https://spliceailookup.broadinstitute.org/#variant="+
					StringUtils.escapeHttp(ucscCtg) + "-" +
					ctx.getStart() +"-"+ctx.getReference().getDisplayString()+"-"+alt.getDisplayString()
					+"&hg="+(isGrch38()?"38":"37")+"&distance=500&mask=1"
					));
			}
		}

	
	
	if(AcidNucleics.isATGC(ctx.getReference()) && (isGrch37() || isGrch38())) {
		for(final Allele alt: ctx.getAlternateAlleles()) {
			if(!AcidNucleics.isATGC(alt)) continue;
			
			//opencravat https://www.opencravat.org/
			urls.add(new LabelledUrlImpl("OpenCravat",
					variantid+"/"+alt.getDisplayString(),
					"https://run.opencravat.org/webapps/variantreport/index.html?chrom="+
					StringUtils.escapeHttp(ucscCtg) + "&pos=" +
					ctx.getStart() +
					"&ref_base=" + ctx.getReference().getDisplayString() +
					"&alt_base=" + alt.getDisplayString()
					+"&assembly="+(isGrch38()?"hg38":"hg19")
					));
			}	
		}
	
	
	//beacon , varsome
	for(int side=0;side<2 && !StringUtils.isBlank(ensemblCtg);++side)
		{
		if(side==0 && !isGrch37()) continue;
		if(side==1 && !isGrch38()) continue;
		for(final Allele alt: ctx.getAlternateAlleles())
			{
			if(ctx.getReference().isSymbolic() || alt.isSymbolic()) continue;
			//https://beacon-network.org/#/search?pos=114267128&chrom=4&allele=A&ref=G&rs=GRCh37
			urls.add(new LabelledUrlImpl("Beacon " + alt.getDisplayString(),
					variantid+"/"+alt.getDisplayString(),
					"https://beacon-network.org/#/search?chrom="+
					StringUtils.escapeHttp(ensemblCtg) +
					"&pos="+ctx.getStart()+
					"&ref="+ ctx.getReference().getDisplayString()+
					"&allele="+ alt.getDisplayString()+
					"&rs="+ (side==0?"GRCh37":"GRCh38")
					));
			urls.add(new LabelledUrlImpl("Varsome " + alt.getDisplayString(),
					variantid+"/"+alt.getDisplayString(),
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
	final String locid = LocatableUtils.toNiceString(loc);
	
	urls.add(new LabelledUrlImpl("IGV",locid,"https://"+ IgvConstants.DEFAULT_HOST +":"+IgvConstants.DEFAULT_PORT + "/goto?locus="+
			StringUtils.escapeHttp(loc.getContig()) + "%3A" +xstart1 +"-"+loc.getEnd()
			));
	
	final String ensemblCtg = toEnsembl.apply(loc.getContig());
	if(isGrch37() && ! StringUtils.isBlank(ensemblCtg)) {
		urls.add(new LabelledUrlImpl("Region Gnomad 2.1",locid,"https://gnomad.broadinstitute.org/region/"+
			StringUtils.escapeHttp(ensemblCtg) + "-" + xstart1 +"-"+ xend1 +"?dataset=gnomad_r2_1"
			));
		urls.add(new LabelledUrlImpl("Region Gnomad SV 2.1",locid,"https://gnomad.broadinstitute.org/region/"+
				StringUtils.escapeHttp(ensemblCtg) + "-" + xstart1 +"-"+ xend1 +"?dataset=gnomad_sv_r2_1"
				));
		urls.add(new LabelledUrlImpl("clinvar 37",locid,"https://www.ncbi.nlm.nih.gov/clinvar/?term="+
				ensemblCtg +
				"%5Bchr%5D+AND+"+ loc.getStart()+"%3A"+ loc.getEnd()+"%5Bchrpos37%5D"
				));
			}
	if(isGrch38() && ! StringUtils.isBlank(ensemblCtg)) {
		urls.add(new LabelledUrlImpl("Region Gnomad 3",locid,"https://gnomad.broadinstitute.org/region/"+
			StringUtils.escapeHttp(ensemblCtg) + "-" + xstart1 +"-"+ xend1 +"?dataset="+GNOMAD_HG38
			));
		
		
		
		urls.add(new LabelledUrlImpl("TogoVar",locid,"https://grch38.togovar.org/?mode=simple&term="+ StringUtils.escapeHttp(ensemblCtg)+":"+loc.getStart()+"-"+loc.getEnd() ));
			
		}
	
	final String ucscCtg =  toUcsc.apply(loc.getContig());
	
	if(isGrch38() && ! StringUtils.isBlank(ucscCtg)) {
		urls.add(new LabelledUrlImpl("Region ResearchAllOfUs",locid,"https://databrowser.researchallofus.org/genomic-variants/"+
			StringUtils.escapeHttp(ucscCtg) + ":" + xstart1 +"-"+ xend1
			));
		urls.add(new LabelledUrlImpl("AF.ukbiobank",locid,"https://afb.ukbiobank.ac.uk/region/"+
				StringUtils.escapeHttp(ucscCtg) + "-" + xstart1 +"-"+ xend1
				));
		}

	
	if(isGrch37() && ! StringUtils.isBlank(ucscCtg)) {		
		if(loc.getLengthOnReference()>1) {
			urls.add(new LabelledUrlImpl("dgv",locid,"http://dgv.tcag.ca/gb2/gbrowse/dgv2_hg19?name="+
					StringUtils.escapeHttp(ucscCtg) + 
					"%3A"+loc.getStart()+"-"+loc.getEnd() + ";search=Search"
					));
			}
		urls.add(new LabelledUrlImpl("pophuman",locid,"https://pophuman.uab.cat/?loc="+
				StringUtils.escapeHttp(ucscCtg) + 
				"%3A"+loc.getStart()+".."+loc.getEnd()
				));
		}
	
	if(isGrch37() && ! StringUtils.isBlank(ensemblCtg)) {
		if(loc.getLengthOnReference()>1) {
			urls.add(new LabelledUrlImpl("decipher",locid,"https://www.deciphergenomics.org/browser#q/grch37:"+
			ensemblCtg + ":" +loc.getStart()+"-"+loc.getEnd()));
			}
		}

	
	
	for(final String build: new String[] {"hg19","hg38","hub_3267197_GCA_009914755.4","mm10","canFam3","canFam4"}) {
		if(StringUtils.isBlank(ucscCtg)) break;
		if(build.equals("hg19") && !SequenceDictionaryUtils.isGRCh37(this.dict)) continue;
		if(build.equals("hg38") && !SequenceDictionaryUtils.isGRCh38(this.dict)) continue;
		if(build.equals("mm10") && !SequenceDictionaryUtils.isGRCm38(this.dict)) continue;
		if(build.equals("canFam3") && !SequenceDictionaryUtils.isCanFam3(this.dict)) continue;
		if(build.equals("canFam4") && !SequenceDictionaryUtils.isCanFam4(this.dict)) continue;
		if(build.equals("hub_3267197_GCA_009914755.4") && !SequenceDictionaryUtils.isCHM13v2(this.dict)) continue;
		
		urls.add(new LabelledUrlImpl("UCSC "+build,
			locid,
			"http://genome.ucsc.edu/cgi-bin/hgTracks?db="+build+"&highlight="+build+"."+
					ucscCtg +
			"%3A"+loc.getStart() +"-"+loc.getEnd() + "&position=" +
			ucscCtg +
			"%3A"+ Math.max(1,loc.getStart()-50) +"-"+(loc.getEnd()+50)
			));
		}

	
	
	if(isGrch38() && ! StringUtils.isBlank(ensemblCtg)) {
		if(loc.getLengthOnReference()>1) {
			urls.add(new LabelledUrlImpl("decipher",locid,"https://www.deciphergenomics.org/browser#q/"+
			ensemblCtg + ":" +loc.getStart()+"-"+loc.getEnd()));
			}
		}

	
	if(isGrch38() && !StringUtils.isBlank(ucscCtg)) {
		urls.add(new LabelledUrlImpl("Decaf Decode Region",
				locid,
				"https://decaf.decode.com/region/"+
				StringUtils.escapeHttp(ucscCtg) +":"+loc.getStart()+"-"+loc.getEnd()
				));
		
		urls.add(new LabelledUrlImpl("VISTA",
				locid,
				"https://enhancer.lbl.gov/vista/browse?filter="+
				StringUtils.escapeHttp(ucscCtg) +"%3A"+loc.getStart()+"-"+loc.getEnd()
				));
		}
	
	
	if(loc.getLengthOnReference()>1) {
		for(int i=0;i< 2;++i) {
			if(i==0 && !isGrch37()) continue;
			if(i==1 && !isGrch38()) continue;
			urls.add(new LabelledUrlImpl("Hi-C",locid,"http://promoter.bx.psu.edu/hi-c/view.php?method=Hi-C&species=human&assembly="+(i==0?"hg19":"hg38")+
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
	else if(loc instanceof Gff3Feature) {
		_gff3(Gff3Feature.class.cast(loc), urls);
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
	private final String target;
	private final String url;
	LabelledUrlImpl(final String label,final String target,final String url) {
		this.label = label;
		this.target = target;
		this.url = url;
		}
	@Override
	public String getLabel() { return label;}
	@Override
	public String getTarget() { return target;}
	@Override
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
	public int compareTo(final LabelledUrl o) {
		int i = this.getDomain().compareTo(o.getDomain());
		if(i!=0) return i;
		i = this.getTarget().compareTo(o.getTarget());
		if(i!=0) return i;
		i = this.getUrl().compareTo(o.getUrl());
		return i;
		}
	
	@Override
	public String toString() {
		return getUrl();
		}
	}
}
