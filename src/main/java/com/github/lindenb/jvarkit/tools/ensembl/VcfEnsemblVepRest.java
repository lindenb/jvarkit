/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.ensembl;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpPost;
import org.apache.http.entity.ByteArrayEntity;
import org.apache.http.entity.ContentType;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;

import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.TeeInputStream;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/**
BEGIN_DOC



## Example

```bash
$ curl -s "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" |\
gunzip -c  | java -jar dist/vcfensemblvep.jar | grep -v '^#' | cut -f 1,2,4,5,8


1	10583	G	A	VEPTRCSQ=processed_transcript|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000456328|A|SO:0001631,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000488147|A|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000541675|A|SO:0001632,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000450305|A|SO:0001631,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000515242|A|SO:0001631,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000538476|A|SO:0001632,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000518655|A|SO:0001631,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000438504|A|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000423562|A|SO:0001632
1	10611	C	G	VEPTRCSQ=processed_transcript|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000456328|G|SO:0001631,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000488147|G|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000541675|G|SO:0001632,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000450305|G|SO:0001631,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000515242|G|SO:0001631,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000538476|G|SO:0001632,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000518655|G|SO:0001631,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000438504|G|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000423562|G|SO:0001632
1	13302	C	T	VEPTRCSQ=processed_transcript|550|550|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000456328|T|SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000488147|T|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000541675|T|SO:0001632,transcribed_unprocessed_pseudogene|342|342|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000450305|T|SO:0001792&SO:0001619,transcribed_unprocessed_pseudogene|543|543|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000515242|T|SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000538476|T|SO:0001632,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000518655|T|SO:0001627&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000438504|T|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000423562|T|SO:0001632
1	13327	G	C	VEPTRCSQ=processed_transcript|575|575|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000456328|C|SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000488147|C|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000541675|C|SO:0001632,transcribed_unprocessed_pseudogene|367|367|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000450305|C|SO:0001792&SO:0001619,transcribed_unprocessed_pseudogene|568|568|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000515242|C|SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000538476|C|SO:0001632,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000518655|C|SO:0001627&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000438504|C|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000423562|C|SO:0001632
1	13957	TC	T	VEPTRCSQ=processed_transcript|1206|1206|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000456328||SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000488147||SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000541675||SO:0001632,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000450305||SO:0001632,transcribed_unprocessed_pseudogene|1199|1199|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000515242||SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000538476||SO:0001632,transcribed_unprocessed_pseudogene|1032|1032|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000518655||SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000438504||SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000423562||SO:0001632
1	13980	T	C	VEPTRCSQ=processed_transcript|1228|1228|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000456328|C|SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000488147|C|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000541675|C|SO:0001632,transcribed_unprocessed_pseudogene|||||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000450305|C|SO:0001632,transcribed_unprocessed_pseudogene|1221|1221|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000515242|C|SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000538476|C|SO:0001632,transcribed_unprocessed_pseudogene|1054|1054|||ENSG00000223972|DDX11L1|HGNC|37102|1|ENST00000518655|C|SO:0001792&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000438504|C|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000423562|C|SO:0001632
1	30923	G	T	VEPTRCSQ=lincRNA|||||ENSG00000243485|MIR1302-10|HGNC|38233|1|ENST00000473358|T|SO:0001627&SO:0001619,lincRNA|||||ENSG00000243485|MIR1302-10|HGNC|38233|1|ENST00000469289|T|SO:0001627&SO:0001619,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000488147|T|SO:0001631,lincRNA|||||ENSG00000237613|FAM138A|HGNC|32334|-1|ENST00000417324|T|SO:0001632,miRNA|||||ENSG00000243485|MIR1302-10|HGNC|38233|1|ENST00000607096|T|SO:0001632,lincRNA|||||ENSG00000237613|FAM138A|HGNC|32334|-1|ENST00000461467|T|SO:0001632,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000538476|T|SO:0001631,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000438504|T|SO:0001631,unprocessed_pseudogene|||||ENSG00000227232|WASH7P|HGNC|38034|-1|ENST00000423562|T|SO:0001631
1	46402	C	CTGT	.
1	47190	G	GA	.
1	51476	T	C	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|C|SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|C|SO:0001631
1	51479	T	A	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|A|SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|A|SO:0001631
1	51914	T	G	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|G|SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|G|SO:0001631
1	51935	C	T	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|T|SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|T|SO:0001631
1	51954	G	C	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|C|SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|C|SO:0001631
1	52058	G	C	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|C|SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|C|SO:0001631
1	52144	T	A	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|A|SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|A|SO:0001631
1	52185	TTAA	T	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647||SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857||SO:0001631
1	52238	T	G	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|G|SO:0001631,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|G|SO:0001631
[INFO/VcfEnsemblVepRest] 2015-03-19 17:19:04 "done: N=20"
[INFO/VcfEnsemblVepRest] 2015-03-19 17:19:04 "Number of Variants:0"
1	53234	CAT	C	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647||SO:0001627&SO:0001619,unprocessed_pseudogene|763|764|||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857||SO:0001792&SO:0001619
1	54353	C	A	VEPTRCSQ=unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000594647|A|SO:0001627&SO:0001619,unprocessed_pseudogene|||||ENSG00000268020|OR4G4P|HGNC|14822|1|ENST00000606857|A|SO:0001632

```

## History

* 2018-02-13: removed XSD, parsing DOM. Added SNPEFF output (but it's incomplete for now...)

END_DOC


 */
@Program(name="vcfensemblvep",
	description="Annotate a VCF with ensembl REST API",
	keywords={"vcf","annotation","rest","ensembl"}
)
public class VcfEnsemblVepRest 
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfEnsemblVepRest.class).make();
	enum OutputFormat {standard,details,base64,snpeff};
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-s","--server"},description="REST server")
	private String server = "http://grch37.rest.ensembl.org";

	@Parameter(names={"-e","--extension"},description="Path extension")
	private String extension = "/vep/homo_sapiens/region";

	@Parameter(names={"-n","--batchSize"},description="batch size. How many variant to send in one HTTP query")
	private int batchSize = 100 ;

	@Parameter(names={"-format","--format"},description="[20180213] Output format")
	private OutputFormat outputFormat = OutputFormat.standard;

	@Parameter(names={"-T","--tee"},description="'Tee' xml response to stderr")
	private boolean teeResponse = false;
	@Parameter(names={"-nofail"},description="[20180213] Do not fail on network error")
	private boolean ignoreNetworkErrors = false;

	
	
	public static final String TAG="VEPTRCSQ";

	private DocumentBuilder documentBuilder;
	private Transformer xmlSerializer;
	private CloseableHttpClient httpClient = null;
	
	
	private static final Set<String> INTERGENIC_CONSEQUENCES_ATTRIBUTES =
			Collections.unmodifiableSet(new LinkedHashSet<String>(Arrays.asList(
			"impact","variant_allele"
			)));
	private static final Set<String> REGULATORY_FEATURES_ATTRIBUTES =
			Collections.unmodifiableSet(new LinkedHashSet<String>(Arrays.asList(
				"biotype","impact","regulatory_feature_id","variant_allele"
			)));
	
	private static final Set<String> COLOCATED_VARIANT_ATTRIBUTES =
			Collections.unmodifiableSet(new LinkedHashSet<String>(Arrays.asList(
			"afr_allele",
			"afr_maf",
			"allele_string",
			"amr_allele",
			"amr_maf",
			"eas_allele",
			"eas_maf",
			"end",
			"eur_allele",
			"eur_maf",
			"gnomad_afr_allele",
			"gnomad_afr_maf",
			"gnomad_allele",
			"gnomad_amr_allele",
			"gnomad_amr_maf",
			"gnomad_asj_allele",
			"gnomad_asj_maf",
			"gnomad_eas_allele",
			"gnomad_eas_maf",
			"gnomad_fin_allele",
			"gnomad_fin_maf",
			"gnomad_maf",
			"gnomad_nfe_allele",
			"gnomad_nfe_maf",
			"gnomad_oth_allele",
			"gnomad_oth_maf",
			"gnomad_sas_allele",
			"gnomad_sas_maf",
			"id",
			"minor_allele",
			"minor_allele_freq",
			"sas_allele",
			"sas_maf",
			"seq_region_name",
			"start",
			"strand"
			)));

	
	private static final Set<String> PREDICTION_ATTRIBUTES=
		Collections.unmodifiableSet(new LinkedHashSet<String>(Arrays.asList(
			"allele_string",
			"assembly_name",
			"end",
			"id",
			"input",
			"most_severe_consequence",
			"seq_region_name",
			"start",
			"strand"
			)));
	
	
	private static final Set<String> TRANSCRIPT_ATTRIBUTES=
		Collections.unmodifiableSet(new LinkedHashSet<String>(Arrays.asList(
			"amino_acids",
			"biotype",
			"canonical",
			"ccds",
			"cdna_end",
			"cdna_start",
			"cds_end",
			"cds_start",
			"codons",
			"distance",
			"exon",
			"gene_id",
			"gene_symbol",
			"gene_symbol_source",
			"hgvsc",
			"hgvsp",
			"hgnc_id",
			"impact",
			"intron",
			"polyphen_prediction",
			"polyphen_score",
			"protein_id",
			"protein_end",
			"protein_start",
			"sift_prediction",
			"sift_score",
			"strand",
			"transcript_id",
			"variant_allele"
			)));
	
	
	private abstract class AbstractAttributeEater
		{
		private final Map<String,String> hash=new HashMap<>();
		protected AbstractAttributeEater(final Set<String> KNOWN_ATTRIBUTES,final Element root) {
			if(root.hasAttributes())
				{
				NamedNodeMap attMap = root.getAttributes();
				for(int i=0;i< attMap.getLength();i++) {
					final Attr att = (Attr)attMap.item(i);
				
					if(!KNOWN_ATTRIBUTES.contains(att.getName()))
						{
						LOG.warn("unknow <"+root.getNodeName()+"> attribute @"+att.getName());
						continue;
						}	
					this.hash.put(att.getName(), att.getValue());
					}
				}
			}
		String get(final String key) {
			return hash.getOrDefault(key, "");
			}
		}
	
	private class IntergenicConsequences extends AbstractAttributeEater
		{
		IntergenicConsequences(final Element root)
			{
			super(INTERGENIC_CONSEQUENCES_ATTRIBUTES,root);
			}
		}
	
	private class RegulatoryFeature extends AbstractAttributeEater
		{
		RegulatoryFeature(final Element root)
			{
			super(REGULATORY_FEATURES_ATTRIBUTES,root);
			}
		
		public String getInfoString()
			{
			switch(VcfEnsemblVepRest.this.outputFormat) {
				case standard:
					return REGULATORY_FEATURES_ATTRIBUTES.stream().
							map(ATT->this.get(ATT)).
							collect(Collectors.joining("|"));
				case details:
					return REGULATORY_FEATURES_ATTRIBUTES.stream().
							map(ATT->ATT+"|"+this.get(ATT)).
							collect(Collectors.joining("|"));
				default: throw new IllegalStateException();
				}
			}
		}

	
	private class ColocatedVariant extends AbstractAttributeEater
		{
		ColocatedVariant(final Element root)
			{
			super(COLOCATED_VARIANT_ATTRIBUTES,root);
			}
		
		public String getInfoString()
			{
			switch(VcfEnsemblVepRest.this.outputFormat) {
				case standard:
					return COLOCATED_VARIANT_ATTRIBUTES.stream().
							map(ATT->this.get(ATT)).
							collect(Collectors.joining("|"));
				case details:
					return COLOCATED_VARIANT_ATTRIBUTES.stream().
							map(ATT->ATT+"|"+this.get(ATT)).
							collect(Collectors.joining("|"));
				default: throw new IllegalStateException();
				}
			}
		}
	
	private class TranscriptConsequence extends AbstractAttributeEater
		{
		private final List<String> consequenceTerms = new ArrayList<>();
		private final List<String> refseq_transcript_ids = new ArrayList<>();
		private final List<String> trembls = new ArrayList<>();
		private final List<String> uniparcs = new ArrayList<>();
		private final List<String> swissprots = new ArrayList<>();
		private final List<String> flags = new ArrayList<>();
		TranscriptConsequence(final Element root) {
			super(TRANSCRIPT_ATTRIBUTES,root);
			
			for(Node c1 = root.getFirstChild();c1!=null;c1=c1.getNextSibling())
				{
				if(c1.getNodeType()!=Node.ELEMENT_NODE) continue;
				final Element e1 = Element.class.cast(c1);
				if(e1.getNodeName().equals("consequence_terms"))
					{
					final String termstr = e1.getTextContent();
					if(StringUtil.isBlank(termstr)) continue;
					this.consequenceTerms.add(termstr);
					}
				else if(e1.getNodeName().equals("refseq_transcript_ids"))
					{
					final String termstr = e1.getTextContent();
					if(StringUtil.isBlank(termstr)) continue;
					this.refseq_transcript_ids.add(termstr);
					}
				else if(e1.getNodeName().equals("trembl"))
					{
					final String termstr = e1.getTextContent();
					if(StringUtil.isBlank(termstr)) continue;
					this.trembls.add(termstr);
					}	
				else if(e1.getNodeName().equals("swissprot"))
					{
					final String termstr = e1.getTextContent();
					if(StringUtil.isBlank(termstr)) continue;
					this.swissprots.add(termstr);
					}
				else if(e1.getNodeName().equals("uniparc"))
					{
					final String termstr = e1.getTextContent();
					if(StringUtil.isBlank(termstr)) continue;
					this.uniparcs.add(termstr);
					}
				else if(e1.getNodeName().equals("flags"))
					{
					final String termstr = e1.getTextContent();
					if(StringUtil.isBlank(termstr)) continue;
					this.flags.add(termstr);
					}
				else
					{
					LOG.warning("unknown element <"+e1.getNodeName()+"> under <"+root.getNodeName()+">");
					}
				}
			}
		
		public String getInfoString()
			{
			switch(VcfEnsemblVepRest.this.outputFormat) {
				case details:
					{
					final StringBuilder sb=new StringBuilder();
					for(final String keyTranscript: TRANSCRIPT_ATTRIBUTES) {
						sb.append(keyTranscript);
						sb.append("|");
						sb.append(this.get(keyTranscript));
						sb.append("|");
						}
					sb.append("so|");
					sb.append(String.join("&", this.consequenceTerms));
					sb.append("|refseq_transcript_ids|");
					sb.append(String.join("&", this.refseq_transcript_ids));
					sb.append("|trembl|");
					sb.append(String.join("&", this.trembls));
					sb.append("|uniparc|");
					sb.append(String.join("&", this.uniparcs));
					sb.append("|swissprot|");
					sb.append(String.join("&", this.swissprots));
					return sb.toString();
					}
				case standard: 
					{
					final StringBuilder sb=new StringBuilder();
					for(final String keyTranscript: TRANSCRIPT_ATTRIBUTES) {
						sb.append(this.get(keyTranscript));
						sb.append("|");
						}
					sb.append(String.join("&", this.consequenceTerms));
					sb.append("|");
					sb.append(String.join("&", this.refseq_transcript_ids));
					return sb.toString();
					}
				case snpeff:
					{
					/* eg:
					 * A|synonymous_variant|LOW|THSD1|ENSG000001820|transcript|ENST0000026132|protein_coding|2/17|c.96G>A|p.Arg32Arg|175/9145|96/3057|32/1018
					 */
					final StringBuilder sb=new StringBuilder();
					sb.append(get("variant_allele"));sb.append("|");
					sb.append(String.join("&", this.consequenceTerms));sb.append("|");
					sb.append("MODIFIER");sb.append("|");/////TODO
					sb.append(get("gene_symbol"));sb.append("|");
					sb.append(get("gene_id"));sb.append("|");
					sb.append(super.hash.containsKey("transcript_id")?"transcript":"");sb.append("|");/////TODO
					sb.append(get("transcript_id"));sb.append("|");
					sb.append(get("biotype"));sb.append("|");
					sb.append("");sb.append("|");//TODO 2/17
					sb.append("");sb.append("|");//TODO c.96G>A
					sb.append("");sb.append("|");//TODO p.Arg32Arg
					sb.append("");sb.append("|");//TODO 175/9145
					sb.append("");sb.append("|");//TODO 96/3057
					sb.append("");//TODO 32/1018
					return sb.toString();
					}
				default: throw new IllegalStateException();
				}
			}

		}
	
	private class EnsVepPrediction extends AbstractAttributeEater
		{
		// input user string
		final List<TranscriptConsequence> transcriptConsequences = new ArrayList<>();
		final List<ColocatedVariant> colocatedVariants = new ArrayList<>();
		final List<RegulatoryFeature> regulatoryFeatures = new ArrayList<>();
		final List<IntergenicConsequences> intergenic_consequences = new ArrayList<>();
		
		EnsVepPrediction(final Element root) {
			super(PREDICTION_ATTRIBUTES,root);
			
			for(Node c1 = root.getFirstChild();c1!=null;c1=c1.getNextSibling())
				{
				if(c1.getNodeType()!=Node.ELEMENT_NODE) continue;
				final Element e1 = Element.class.cast(c1);
				if(e1.getNodeName().equals("transcript_consequences"))
					{
					this.transcriptConsequences.add(parseTranscriptConsequence(e1));
					}
				else if(e1.getNodeName().equals("colocated_variants"))
					{
					this.colocatedVariants.add(parseColocatedVariant(e1));
					}
				else if(e1.getNodeName().equals("regulatory_feature_consequences"))
					{
					this.regulatoryFeatures.add(new RegulatoryFeature(e1));
					}
				else if(e1.getNodeName().equals("intergenic_consequences"))
					{
					this.intergenic_consequences.add(new IntergenicConsequences(e1));
					}
				else
					{
					LOG.warn("unknow element <"+e1.getNodeName()+">  under <"+root.getNodeName()+">");
					}
				}
			}
		
		String getInput() { return this.get("input");}
		
		public List<String> getInfoStringList()
			{
			switch(VcfEnsemblVepRest.this.outputFormat) {
				case standard: 
				case snpeff:
				case details:
					{
					return transcriptConsequences.stream().
							map(T->T.getInfoString()).
							collect(Collectors.toList());
					}
				default: throw new IllegalStateException();
				}
			}
		public List<String> getCollocatedVariantStringList()
			{
			switch(VcfEnsemblVepRest.this.outputFormat) {
				case standard: 
				case details:
					{
					return this.colocatedVariants.stream().
							map(T->T.getInfoString()).
							collect(Collectors.toList());
					}
				case snpeff:
					{
					return Collections.emptyList();
					}
				default: throw new IllegalStateException();
				}
			}
		public List<String> getRegulationStringList()
			{
			switch(VcfEnsemblVepRest.this.outputFormat) {
				case standard: 
				case details:
					{
					return this.regulatoryFeatures.stream().
							map(T->T.getInfoString()).
							collect(Collectors.toList());
					}
				case snpeff:
					{
					return Collections.emptyList();
					}
				default: throw new IllegalStateException();
				}
			}

		}
	
	private ColocatedVariant parseColocatedVariant(final Element root) {
		return new ColocatedVariant(root);
		}
	
	private TranscriptConsequence parseTranscriptConsequence(final Element root) {
		final TranscriptConsequence pred = new TranscriptConsequence(root);
		
		return pred;
		}
	
	private EnsVepPrediction parseVepPrediction(final Element root) {
		final EnsVepPrediction pred = new EnsVepPrediction(root);
		
		return pred;
		}
	
	private List<EnsVepPrediction> parseVepPredictions(final Document dom) {
		final List<EnsVepPrediction> preds = new ArrayList<>();
		final Element root= dom.getDocumentElement();
		if(root==null) {
			if(ignoreNetworkErrors) {
				LOG.error("empty dom");
				return Collections.emptyList();
				}
			throw new IllegalStateException("empty dom");
			}
		if(!root.getNodeName().equals("opt")) {
			final String msg = "root is not <opt> but <"+root.getNodeName()+">";
			if(ignoreNetworkErrors) {
				LOG.error( msg);
				return Collections.emptyList();
				}
			throw new IllegalStateException( msg);
		}
		for(Node c1 = root.getFirstChild();c1!=null;c1=c1.getNextSibling())
			{
			if(c1.getNodeType()!=Node.ELEMENT_NODE) continue;
			final Element e1 = Element.class.cast(c1);
			if(e1.getNodeName().equals("data"))
				{
				preds.add(parseVepPrediction(e1));
				}	
			}
		return preds;
		}
	
	private static String createInputContext(final VariantContext ctx)
		{
		final StringBuilder sb=new StringBuilder();
		sb.append(ctx.getContig()).
			append(" ").
			append(ctx.getStart()).
			append(" ").
			append(!ctx.hasID()?".":ctx.getID()).
			append(" ").
			append(ctx.getReference().getBaseString()).
			append(" ")
			;
		final List<Allele> alts=ctx.getAlternateAlleles();
		if(alts.isEmpty())
			{
			sb.append(".");
			}
		else
			{
			for(int j=0;j< alts.size();++j )
				{
				if(j>0) sb.append(",");
				sb.append(alts.get(j).getBaseString());
				}
			}
		sb.append(" . . .");
		return sb.toString();
		}
		
	private long lastMillisec=-1L;
	
	/** send a pool of variants to VEP, returns the DOM document */
	private Document callVepToDom(final List<VariantContext> contexts,boolean xml_answer) throws IOException
		{
		LOG.info("Running VEP "+contexts.size());
		InputStream response =null;
		HttpPost httpPost = null;
		try {
		    if ( this.lastMillisec!=-1L && this.lastMillisec+ 5000<  System.currentTimeMillis())
		    	{
		    	LOG.debug("waiting");
		    	try {Thread.sleep(1000);} catch(Exception err){}
		    	}
				 
		    httpPost = new HttpPost(this.server + this.extension);
			 
			 final StringBuilder queryb=new StringBuilder();
			 queryb.append("{ \"variants\" : [");
			 for(int i=0;i< contexts.size();++i)
			 	{
				final VariantContext ctx=contexts.get(i);
				if(i>0) queryb.append(",");
				queryb.append("\"").
					append(createInputContext(ctx)).
					append("\"");
			 	}
			 queryb.append("]");
			 for(final String s: new String[]{
					 "canonical","ccds","domains","hgvs","numbers",
					 "protein","xref_refseq","tsl","uniprot"})
			 	{
				 queryb.append(",\"").append(s).append("\":1");
			 	}
			 queryb.append("}");
			 final byte postBody[] = queryb.toString().getBytes();

			 httpPost.setHeader("Content-Type",ContentType.APPLICATION_JSON.getMimeType());
			 httpPost.setHeader("Accept",ContentType.TEXT_XML.getMimeType());
			 //httpPost.setHeader("Content-Length", Integer.toString(postBody.length));
			 httpPost.setEntity(new ByteArrayEntity(postBody, ContentType.APPLICATION_JSON));
			 
			 
			 
			 final CloseableHttpResponse httpResponse = httpClient.execute(httpPost);
			 
			 int responseCode = httpResponse.getStatusLine().getStatusCode();
			  
			 if(responseCode != 200)
			 	{
				throw new RuntimeIOException("Response code was not 200. Detected response was "+responseCode);
			 	}

			 
			 //response = new TeeInputStream( httpConnection.getInputStream(),System.err,false);
			 response =httpResponse.getEntity().getContent();
			 if(this.teeResponse)
				 {
				 stderr().println(queryb);
				 response = new TeeInputStream(response,stderr(),false);
				 }
			 
			final Document dom = documentBuilder.parse(response);
			return dom;
			} 
		catch (final Throwable err)
			{
			if(this.ignoreNetworkErrors) {
				LOG.error(err);
				return documentBuilder.newDocument();
				}
			throw new IOException(err);
			}
		finally
			{
			CloserUtil.close(response);
			if(httpPost!=null) httpPost.releaseConnection();
			this.lastMillisec = System.currentTimeMillis(); 
			}
		}
	
	private List<EnsVepPrediction> vep(final List<VariantContext> contexts) throws IOException
		{
		return parseVepPredictions(callVepToDom(contexts,false));
		}
	private Document vepxml(List<VariantContext> contexts) throws IOException
		{
		return callVepToDom(contexts,true);
		}
	
	@Override
	public int doWork(final List<String> args) {
	try {
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		this.documentBuilder=dbf.newDocumentBuilder();		
		TransformerFactory trf=TransformerFactory.newInstance();
		this.xmlSerializer = trf.newTransformer();
		
		/** create http client */
		this.httpClient = HttpClients.createSystem();//createDefault();
		
		return doVcfToVcf(args, this.outputFile);
		}
	catch(Exception err) {
		LOG.error(err);
		return -1;
		}
	finally 
		{
		CloserUtil.close(this.httpClient);
		this.httpClient=null;
		}
	}
		
	
	@Override
	protected int doVcfToVcf(final String inputName, final VcfIterator vcfIn, final VariantContextWriter out) {
	    try {
		final java.util.Base64.Encoder  base64Encoder;
		final VCFHeader header=vcfIn.getHeader();
		final List<VariantContext> buffer=new ArrayList<>(this.batchSize+1);
		final VCFHeader h2= new VCFHeader(header);
		addMetaData(h2);
		final String tagName1;
		final String tagName2;
		final String tagReg;
		switch(this.outputFormat)
			{
			case base64:
				{
				tagName1 = TAG+"_base64";
				tagName2 = null;
				tagReg = null;
				base64Encoder = java.util.Base64.getEncoder();
				h2.addMetaDataLine(new VCFInfoHeaderLine(
						tagName1,
						1,
						VCFHeaderLineType.String,
						"VEP xml answer encoded as base 64"
						));
				break;
				}
			case snpeff:
				{
				tagName1 = "ANN";
				tagName2 = null;
				tagReg = null;
				base64Encoder = null;
				h2.addMetaDataLine(new VCFInfoHeaderLine(
						tagName1,
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'"
						));
				LOG.warn("ANN/snpeff output is not fully implemented");
				break;
				}
			case details:
				{
				tagName1 = TAG+"X";
				tagName2 = TAG+ "X_COLOCATED_VARIANTS";
				tagReg = TAG+ "X_REG";
				base64Encoder = null;
				h2.addMetaDataLine(new VCFInfoHeaderLine(
						tagName1,
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"VEP Transcript Consequences. Format :("+TRANSCRIPT_ATTRIBUTES.stream().map(S->S.toUpperCase()+"|"+S.toLowerCase()).collect(Collectors.joining("|"))+"|SO|so|REFSEQ|refseq)"
						));
				h2.addMetaDataLine(new VCFInfoHeaderLine(
						tagName2,
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"VEP Colocated variants.  Format :("+COLOCATED_VARIANT_ATTRIBUTES.stream().map(S->S.toUpperCase()+"|"+S.toLowerCase()).collect(Collectors.joining("|"))+")"
						));
				h2.addMetaDataLine(new VCFInfoHeaderLine(
						tagReg,
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"VEP Regulatory Features  Format :("+REGULATORY_FEATURES_ATTRIBUTES.stream().map(S->S.toUpperCase()+"|"+S.toLowerCase()).collect(Collectors.joining("|"))+")"
						));
				break;
				}
			case standard:
				{
				tagName1 = TAG;
				tagName2 = TAG+ "_COLOCATED_VARIANTS";
				tagReg = TAG+"_REG";
				base64Encoder = null;
				h2.addMetaDataLine(new VCFInfoHeaderLine(
						tagName1,
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"VEP Transcript Consequences. Format :("+String.join("|",TRANSCRIPT_ATTRIBUTES)+"|so_acns)"
						));
				h2.addMetaDataLine(new VCFInfoHeaderLine(
						tagName2,
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"VEP Colocated variants.  Format :("+String.join("|",COLOCATED_VARIANT_ATTRIBUTES)+")"
						));
				h2.addMetaDataLine(new VCFInfoHeaderLine(
						tagReg,
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"VEP Colocated variants.  Format :("+String.join("|",REGULATORY_FEATURES_ATTRIBUTES)+")"
						));
				break;
				}
			default: throw new IllegalArgumentException();
			}
		
		out.writeHeader(h2);
		final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header).logger(LOG);
		for(;;)
			{
			VariantContext ctx=null;
			if(vcfIn.hasNext())
				{
				buffer.add((ctx=progress.watch(vcfIn.next())));
				}
			if(ctx==null || buffer.size()>=this.batchSize)
				{
				if(!buffer.isEmpty())
					{
					switch(this.outputFormat)
						{
						case base64:
							{
							final Document opt = vepxml(buffer);
							final Element root= opt.getDocumentElement();
							if(!root.getNodeName().equals("opt"))
								throw new IOException("Bad root node "+root.getNodeName());
							
							for(final VariantContext ctx2:buffer)
								{
								String inputStr = createInputContext(ctx2);							
								Document newdom=null;
								
								//loop over <data/>
								for(Node dataNode =root.getFirstChild();
										dataNode!=null;
										dataNode=dataNode.getNextSibling())
									{
									if(dataNode.getNodeType()!=Node.ELEMENT_NODE) continue;
									final Attr att = Element.class.cast(dataNode).getAttributeNode("input");
									if(att==null)
										{
										LOG.warn("no @input in <data/>");
										continue;
										}
	
									if(!att.getValue().equals(inputStr)) continue;
									if(newdom==null)
										{
										newdom = this.documentBuilder.newDocument();
										newdom.appendChild(newdom.createElement("opt"));
										}
									newdom.getDocumentElement().appendChild(newdom.importNode(dataNode, true));
									}
								if(newdom==null)
									{
									LOG.warn("No Annotation found for "+inputStr);
									out.add(ctx2);
									continue;
									}
								final StringWriter sw=new StringWriter();
								try {
									this.xmlSerializer.transform(
											new DOMSource(newdom),
											new StreamResult(sw)
											);
									} 
								catch(final TransformerException err)
									{
									throw new IOException(err);
									}
								final VariantContextBuilder vcb=new VariantContextBuilder(ctx2);
								vcb.attribute(TAG,base64Encoder.encodeToString(sw.toString().getBytes()).
										replaceAll("[\\s=]", ""));
								out.add(vcb.make());
								}
							break;
							}
						case snpeff:
						case details:
						default:
							{
							final List<EnsVepPrediction> predlist = vep(buffer);
							for(final VariantContext ctx2:buffer)
								{
								final String inputStr = createInputContext(ctx2);
								final Optional<EnsVepPrediction> optmydata = 
										predlist.stream().
										filter(P->inputStr.equals(P.getInput())).
										findFirst()
										;
								if(!optmydata.isPresent())
									{
									LOG.info("No Annotation found for "+inputStr);
									out.add(ctx2);
									continue;
									}
								final List<String> infoList = optmydata.get().getInfoStringList();
								final List<String> otherVariantsList=  optmydata.get().getCollocatedVariantStringList();
								final List<String> regulatoryList=  optmydata.get().getRegulationStringList();
								
								final VariantContextBuilder vcb=new VariantContextBuilder(ctx2);
								
								if(!infoList.isEmpty()) {
									vcb.attribute(tagName1, infoList);
									}
								if(!otherVariantsList.isEmpty()) {
									vcb.attribute(tagName2, otherVariantsList);
									}
								if(!regulatoryList.isEmpty()) {
									vcb.attribute(tagReg, regulatoryList);
									}
								out.add(vcb.make());
								}
							break;
							}//end of not(XML base 64)
						}//end if switch output format
					} // end of if buffer is not empty
				if(ctx==null) break;
				buffer.clear();
				}
			}
		progress.finish();
		return RETURN_OK;
	    } 
	catch(final Exception err)
    	{
    	LOG.error(err);
    	return -1;
    	}

	}
	
	
	public static void main(String[] args) {
		new VcfEnsemblVepRest().instanceMainWithExit(args);
	}
	}
