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
package com.github.lindenb.jvarkit.tools.ensembl;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
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
import javax.xml.transform.Source;
import javax.xml.transform.Templates;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.stream.StreamSource;

import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpPost;
import org.apache.http.entity.ByteArrayEntity;
import org.apache.http.entity.ContentType;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.TeeInputStream;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC



## Example

```bash

 $ wget -q -O - "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"  |\
  	gunzip -c  |\
  	java -jar dist/vcfensemblvep.jar  -n 10 | grep -v '^#' | cut -f 1,2,4,5,8 | head

1	10583	G	A	AA=.;AC=314;AF=0.14;AFR_AF=0.04;AMR_AF=0.17;AN=2184;ASN_AF=0.13;AVGPOST=0.7707;ERATE=0.0161;EUR_AF=0.21;LDAF=0.2327;RSQ=0.4319;SNPSOURCE=LOWCOV;THETA=0.0046;VEPREST_COLOCATED_VARIANTS=allele_string|G/A|end|10583|id|rs58108140|seq_region_name|1|start|10583|strand|1;VEPREST_TRANSCRIPT=biotype|unprocessed_pseudogene|distance|3780|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000423562|variant_allele|A|consequences|downstream_gene_variant,biotype|unprocessed_pseudogene|canonical|1|distance|3780|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000438504|variant_allele|A|consequences|downstream_gene_variant,biotype|transcribed_unprocessed_pseudogene|distance|1427|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|impact|MODIFIER|strand|1|transcript_id|ENST00000450305|variant_allele|A|consequences|upstream_gene_variant,biotype|processed_transcript|canonical|1|distance|1286|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|impact|MODIFIER|strand|1|transcript_id|ENST00000456328|variant_allele|A|consequences|upstream_gene_variant,biotype|unprocessed_pseudogene|distance|3821|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000488147|variant_allele|A|consequences|downstream_gene_variant,biotype|transcribed_unprocessed_pseudogene|distance|1289|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|impact|MODIFIER|strand|1|transcript_id|ENST00000515242|variant_allele|A|consequences|upstream_gene_variant,biotype|transcribed_unprocessed_pseudogene|distance|1291|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|impact|MODIFIER|strand|1|transcript_id|ENST00000518655|variant_allele|A|consequences|upstream_gene_variant,biotype|unprocessed_pseudogene|distance|3828|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000538476|variant_allele|A|consequences|downstream_gene_variant,biotype|unprocessed_pseudogene|distance|3780|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000541675|variant_allele|A|consequences|downstream_gene_variant;VT=SNP
1	10611	C	G	AA=.;AC=41;AF=0.02;AFR_AF=0.01;AMR_AF=0.03;AN=2184;ASN_AF=0.01;AVGPOST=0.9330;ERATE=0.0048;EUR_AF=0.02;LDAF=0.0479;RSQ=0.3475;SNPSOURCE=LOWCOV;THETA=0.0077;VEPREST_COLOCATED_VARIANTS=allele_string|C/G|end|10611|id|rs189107123|seq_region_name|1|start|10611|strand|1;VEPREST_TRANSCRIPT=biotype|unprocessed_pseudogene|distance|3752|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000423562|variant_allele|G|consequences|downstream_gene_variant,biotype|unprocessed_pseudogene|canonical|1|distance|3752|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000438504|variant_allele|G|consequences|downstream_gene_variant,biotype|transcribed_unprocessed_pseudogene|distance|1399|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|impact|MODIFIER|strand|1|transcript_id|ENST00000450305|variant_allele|G|consequences|upstream_gene_variant,biotype|processed_transcript|canonical|1|distance|1258|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|impact|MODIFIER|strand|1|transcript_id|ENST00000456328|variant_allele|G|consequences|upstream_gene_variant,biotype|unprocessed_pseudogene|distance|3793|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000488147|variant_allele|G|consequences|downstream_gene_variant,biotype|transcribed_unprocessed_pseudogene|distance|1261|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|impact|MODIFIER|strand|1|transcript_id|ENST00000515242|variant_allele|G|consequences|upstream_gene_variant,biotype|transcribed_unprocessed_pseudogene|distance|1263|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|impact|MODIFIER|strand|1|transcript_id|ENST00000518655|variant_allele|G|consequences|upstream_gene_variant,biotype|unprocessed_pseudogene|distance|3800|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000538476|variant_allele|G|consequences|downstream_gene_variant,biotype|unprocessed_pseudogene|distance|3752|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000541675|variant_allele|G|consequences|downstream_gene_variant;VT=SNP
1	13302	C	T	AA=.;AC=249;AF=0.11;AFR_AF=0.21;AMR_AF=0.08;AN=2184;ASN_AF=0.02;AVGPOST=0.8895;ERATE=0.0058;EUR_AF=0.14;LDAF=0.1573;RSQ=0.6281;SNPSOURCE=LOWCOV;THETA=0.0048;VEPREST_COLOCATED_VARIANTS=allele_string|C/G/T|end|13302|id|rs75241669|seq_region_name|1|start|13302|strand|1;VEPREST_TRANSCRIPT=biotype|unprocessed_pseudogene|distance|1061|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000423562|variant_allele|T|consequences|downstream_gene_variant,biotype|unprocessed_pseudogene|canonical|1|distance|1061|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000438504|variant_allele|T|consequences|downstream_gene_variant,biotype|transcribed_unprocessed_pseudogene|cdna_end|342|cdna_start|342|exon|5/6|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|hgvsc|ENST00000450305.2:n.342C>T|impact|MODIFIER|strand|1|transcript_id|ENST00000450305|variant_allele|T|consequences|non_coding_transcript_exon_variant,biotype|processed_transcript|canonical|1|cdna_end|550|cdna_start|550|exon|3/3|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|hgvsc|ENST00000456328.2:n.550C>T|impact|MODIFIER|strand|1|transcript_id|ENST00000456328|variant_allele|T|consequences|non_coding_transcript_exon_variant,biotype|unprocessed_pseudogene|distance|1102|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000488147|variant_allele|T|consequences|downstream_gene_variant,biotype|transcribed_unprocessed_pseudogene|cdna_end|543|cdna_start|543|exon|3/3|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|hgvsc|ENST00000515242.2:n.543C>T|impact|MODIFIER|strand|1|transcript_id|ENST00000515242|variant_allele|T|consequences|non_coding_transcript_exon_variant,biotype|transcribed_unprocessed_pseudogene|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|hgvsc|ENST00000518655.2:n.482-101C>T|impact|MODIFIER|intron|2/3|strand|1|transcript_id|ENST00000518655|variant_allele|T|consequences|intron_variant&non_coding_transcript_variant,biotype|unprocessed_pseudogene|distance|1109|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000538476|variant_allele|T|consequences|downstream_gene_variant,biotype|unprocessed_pseudogene|distance|1061|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000541675|variant_allele|T|consequences|downstream_gene_variant;VT=SNP
1	13327	G	C	AA=.;AC=59;AF=0.03;AFR_AF=0.02;AMR_AF=0.03;AN=2184;ASN_AF=0.02;AVGPOST=0.9698;ERATE=0.0012;EUR_AF=0.04;LDAF=0.0359;RSQ=0.6482;SNPSOURCE=LOWCOV;THETA=0.0204;VEPREST_TRANSCRIPT=biotype|unprocessed_pseudogene|distance|1036|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000423562|variant_allele|C|consequences|downstream_gene_variant,biotype|unprocessed_pseudogene|canonical|1|distance|1036|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000438504|variant_allele|C|consequences|downstream_gene_variant,biotype|transcribed_unprocessed_pseudogene|cdna_end|367|cdna_start|367|exon|5/6|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|hgvsc|ENST00000450305.2:n.367G>C|impact|MODIFIER|strand|1|transcript_id|ENST00000450305|variant_allele|C|consequences|non_coding_transcript_exon_variant,biotype|processed_transcript|canonical|1|cdna_end|575|cdna_start|575|exon|3/3|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|hgvsc|ENST00000456328.2:n.575G>C|impact|MODIFIER|strand|1|transcript_id|ENST00000456328|variant_allele|C|consequences|non_coding_transcript_exon_variant,biotype|unprocessed_pseudogene|distance|1077|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000488147|variant_allele|C|consequences|downstream_gene_variant,biotype|transcribed_unprocessed_pseudogene|cdna_end|568|cdna_start|568|exon|3/3|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|hgvsc|ENST00000515242.2:n.568G>C|impact|MODIFIER|strand|1|transcript_id|ENST00000515242|variant_allele|C|consequences|non_coding_transcript_exon_variant,biotype|transcribed_unprocessed_pseudogene|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|hgvsc|ENST00000518655.2:n.482-76G>C|impact|MODIFIER|intron|2/3|strand|1|transcript_id|ENST00000518655|variant_allele|C|consequences|intron_variant&non_coding_transcript_variant,biotype|unprocessed_pseudogene|distance|1084|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000538476|variant_allele|C|consequences|downstream_gene_variant,biotype|unprocessed_pseudogene|distance|1036|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000541675|variant_allele|C|consequences|downstream_gene_variant;VT=SNP
1	13957	TC	T	AA=TC;AC=35;AF=0.02;AFR_AF=0.02;AMR_AF=0.02;AN=2184;ASN_AF=0.01;AVGPOST=0.8711;ERATE=0.0065;EUR_AF=0.02;LDAF=0.0788;RSQ=0.2501;THETA=0.0100;VEPREST_COLOCATED_VARIANTS=allele_string|C/-|end|13958|id|rs201747181|seq_region_name|1|start|13958|strand|1;VEPREST_TRANSCRIPT=biotype|unprocessed_pseudogene|distance|405|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000423562|variant_allele|-|consequences|downstream_gene_variant,biotype|unprocessed_pseudogene|canonical|1|distance|405|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000438504|variant_allele|-|consequences|downstream_gene_variant,biotype|transcribed_unprocessed_pseudogene|distance|288|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|impact|MODIFIER|strand|1|transcript_id|ENST00000450305|variant_allele|-|consequences|downstream_gene_variant,biotype|processed_transcript|canonical|1|cdna_end|1206|cdna_start|1206|exon|3/3|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|hgvs_offset|4|hgvsc|ENST00000456328.2:n.1210del|impact|MODIFIER|strand|1|transcript_id|ENST00000456328|variant_allele|-|consequences|non_coding_transcript_exon_variant,biotype|unprocessed_pseudogene|distance|446|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000488147|variant_allele|-|consequences|downstream_gene_variant,biotype|transcribed_unprocessed_pseudogene|cdna_end|1199|cdna_start|1199|exon|3/3|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|hgvs_offset|4|hgvsc|ENST00000515242.2:n.1203del|impact|MODIFIER|strand|1|transcript_id|ENST00000515242|variant_allele|-|consequences|non_coding_transcript_exon_variant,biotype|transcribed_unprocessed_pseudogene|cdna_end|1032|cdna_start|1032|exon|4/4|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|hgvs_offset|4|hgvsc|ENST00000518655.2:n.1036del|impact|MODIFIER|strand|1|transcript_id|ENST00000518655|variant_allele|-|consequences|non_coding_transcript_exon_variant,biotype|unprocessed_pseudogene|distance|453|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000538476|variant_allele|-|consequences|downstream_gene_variant,biotype|unprocessed_pseudogene|distance|405|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000541675|variant_allele|-|consequences|downstream_gene_variant;VT=INDEL
1	13980	T	C	AA=.;AC=45;AF=0.02;AFR_AF=0.01;AMR_AF=0.02;AN=2184;ASN_AF=0.02;AVGPOST=0.9221;ERATE=0.0034;EUR_AF=0.02;LDAF=0.0525;RSQ=0.3603;SNPSOURCE=LOWCOV;THETA=0.0139;VEPREST_TRANSCRIPT=biotype|unprocessed_pseudogene|distance|383|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000423562|variant_allele|C|consequences|downstream_gene_variant,biotype|unprocessed_pseudogene|canonical|1|distance|383|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000438504|variant_allele|C|consequences|downstream_gene_variant,biotype|transcribed_unprocessed_pseudogene|distance|310|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|impact|MODIFIER|strand|1|transcript_id|ENST00000450305|variant_allele|C|consequences|downstream_gene_variant,biotype|processed_transcript|canonical|1|cdna_end|1228|cdna_start|1228|exon|3/3|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|hgvsc|ENST00000456328.2:n.1228T>C|impact|MODIFIER|strand|1|transcript_id|ENST00000456328|variant_allele|C|consequences|non_coding_transcript_exon_variant,biotype|unprocessed_pseudogene|distance|424|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000488147|variant_allele|C|consequences|downstream_gene_variant,biotype|transcribed_unprocessed_pseudogene|cdna_end|1221|cdna_start|1221|exon|3/3|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|hgvsc|ENST00000515242.2:n.1221T>C|impact|MODIFIER|strand|1|transcript_id|ENST00000515242|variant_allele|C|consequences|non_coding_transcript_exon_variant,biotype|transcribed_unprocessed_pseudogene|cdna_end|1054|cdna_start|1054|exon|4/4|gene_id|ENSG00000223972|gene_symbol|DDX11L1|gene_symbol_source|HGNC|hgnc_id|37102|hgvsc|ENST00000518655.2:n.1054T>C|impact|MODIFIER|strand|1|transcript_id|ENST00000518655|variant_allele|C|consequences|non_coding_transcript_exon_variant,biotype|unprocessed_pseudogene|distance|431|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000538476|variant_allele|C|consequences|downstream_gene_variant,biotype|unprocessed_pseudogene|distance|383|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000541675|variant_allele|C|consequences|downstream_gene_variant;VT=SNP
1	30923	G	T	AA=T;AC=1584;AF=0.73;AFR_AF=0.48;AMR_AF=0.80;AN=2184;ASN_AF=0.89;AVGPOST=0.7335;ERATE=0.0183;EUR_AF=0.73;LDAF=0.6576;RSQ=0.5481;SNPSOURCE=LOWCOV;THETA=0.0162;VEPREST_COLOCATED_VARIANTS=afr_allele|T|afr_maf|0.6687|allele_string|G/T|amr_allele|T|amr_maf|0.9164|eas_allele|T|eas_maf|0.996|end|30923|eur_allele|T|eur_maf|0.9364|id|rs806731|minor_allele|G|minor_allele_freq|0.1276|sas_allele|T|sas_maf|0.9233|seq_region_name|1|start|30923|strand|1;VEPREST_TRANSCRIPT=biotype|lincRNA|canonical|1|distance|3631|gene_id|ENSG00000237613|gene_symbol|FAM138A|gene_symbol_source|HGNC|hgnc_id|32334|impact|MODIFIER|strand|-1|transcript_id|ENST00000417324|variant_allele|T|consequences|downstream_gene_variant,biotype|unprocessed_pseudogene|distance|1553|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000423562|variant_allele|T|consequences|upstream_gene_variant,biotype|unprocessed_pseudogene|canonical|1|distance|1553|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000438504|variant_allele|T|consequences|upstream_gene_variant,biotype|lincRNA|distance|4322|gene_id|ENSG00000237613|gene_symbol|FAM138A|gene_symbol_source|HGNC|hgnc_id|32334|impact|MODIFIER|strand|-1|transcript_id|ENST00000461467|variant_allele|T|consequences|downstream_gene_variant,biotype|lincRNA|gene_id|ENSG00000243485|gene_symbol|MIR1302-10|gene_symbol_source|HGNC|hgnc_id|38233|hgvsc|ENST00000469289.1:n.402-53G>T|impact|MODIFIER|intron|1/1|strand|1|transcript_id|ENST00000469289|variant_allele|T|consequences|intron_variant&non_coding_transcript_variant,biotype|lincRNA|canonical|1|gene_id|ENSG00000243485|gene_symbol|MIR1302-10|gene_symbol_source|HGNC|hgnc_id|38233|hgvsc|ENST00000473358.1:n.591-53G>T|impact|MODIFIER|intron|2/2|strand|1|transcript_id|ENST00000473358|variant_allele|T|consequences|intron_variant&non_coding_transcript_variant,biotype|unprocessed_pseudogene|distance|1353|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000488147|variant_allele|T|consequences|upstream_gene_variant,biotype|unprocessed_pseudogene|distance|1117|gene_id|ENSG00000227232|gene_symbol|WASH7P|gene_symbol_source|HGNC|hgnc_id|38034|impact|MODIFIER|strand|-1|transcript_id|ENST00000538476|variant_allele|T|consequences|upstream_gene_variant,biotype|miRNA|distance|420|gene_id|ENSG00000243485|gene_symbol|MIR1302-10|gene_symbol_source|HGNC|hgnc_id|38233|impact|MODIFIER|strand|1|transcript_id|ENST00000607096|variant_allele|T|consequences|downstream_gene_variant;VT=SNP
1	46402	C	CTGT	AA=.;AC=8;AF=0.0037;AFR_AF=0.01;AN=2184;ASN_AF=0.0017;AVGPOST=0.8325;ERATE=0.0072;LDAF=0.0903;RSQ=0.0960;THETA=0.0121;VEPREST_COLOCATED_VARIANTS=allele_string|-/TGT|end|46402|id|rs199681827|seq_region_name|1|start|46403|strand|1;VEPREST_INTERGENIC_CSQ=impact|MODIFIER|variant_allele|TGT|consequences|intergenic_variant;VT=INDEL
1	47190	G	GA	AA=G;AC=29;AF=0.01;AFR_AF=0.06;AMR_AF=0.0028;AN=2184;AVGPOST=0.9041;ERATE=0.0041;LDAF=0.0628;RSQ=0.2883;THETA=0.0153;VEPREST_COLOCATED_VARIANTS=allele_string|-/A|end|47190|id|rs200430748|seq_region_name|1|start|47191|strand|1;VEPREST_INTERGENIC_CSQ=impact|MODIFIER|variant_allele|A|consequences|intergenic_variant;VT=INDEL
1	51476	T	C	AA=C;AC=18;AF=0.01;AFR_AF=0.01;AMR_AF=0.01;AN=2184;ASN_AF=0.01;AVGPOST=0.9819;ERATE=0.0021;EUR_AF=0.01;LDAF=0.0157;RSQ=0.5258;SNPSOURCE=LOWCOV;THETA=0.0103;VEPREST_TRANSCRIPT=biotype|unprocessed_pseudogene|distance|1573|gene_id|ENSG00000268020|gene_symbol|OR4G4P|gene_symbol_source|HGNC|hgnc_id|14822|impact|MODIFIER|strand|1|transcript_id|ENST00000594647|variant_allele|C|consequences|upstream_gene_variant,biotype|unprocessed_pseudogene|canonical|1|distance|997|gene_id|ENSG00000268020|gene_symbol|OR4G4P|gene_symbol_source|HGNC|hgnc_id|14822|impact|MODIFIER|strand|1|transcript_id|ENST00000606857|variant_allele|C|consequences|upstream_gene_variant;VT=SNP

```

## History

* 2018-02-13: removed XSD, parsing DOM. Added SNPEFF output (but it's incomplete for now...)

END_DOC


 */
@Program(name="vcfensemblvep",
	description="Annotate a VCF with ensembl REST API",
	keywords={"vcf","annotation","rest","ensembl","xml","xslt","xsl"}
)
public class VcfEnsemblVepRest 
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfEnsemblVepRest.class).make();
	private static final String DEFAULT_XSLT_TEMPLATE = "<?xml version=\'1.0\' encoding=\"UTF-8\"?>\n" + 
			"<xsl:stylesheet xmlns:xsl=\'http://www.w3.org/1999/XSL/Transform\' version=\'1.0\'>\n" + 
			"<xsl:output method=\"text\"/>\n" + 
			"\n" + 
			"<xsl:template match=\"/\">\n" + 
			"  <xsl:apply-templates select=\"*\"/>\n" + 
			"</xsl:template>\n" + 
			"\n" + 
			"<xsl:template match=\"opt\">\n" + 
			"  <xsl:apply-templates select=\"data\"/>\n" + 
			"</xsl:template>\n" + 
			"\n" + 
			"<xsl:template match=\"data\">\n" + 
			" <xsl:if test=\"transcript_consequences\">\n" + 
			" <xsl:text>VEPREST_TRANSCRIPT=</xsl:text>\n" + 
			" <xsl:for-each select=\"transcript_consequences\">\n" + 
			"   <xsl:if test=\"position()&gt;1\">,</xsl:if>\n" + 
			"   <xsl:apply-templates select=\".\"/>\n" + 
			" </xsl:for-each>\n" + 
			"<xsl:text>\n" + 
			"</xsl:text>\n" + 
			" </xsl:if>\n" + 
			" \n" + 
			" <xsl:if test=\"colocated_variants\">\n" + 
			" <xsl:text>VEPREST_COLOCATED_VARIANTS=</xsl:text>\n" + 
			" <xsl:for-each select=\"colocated_variants\">\n" + 
			"   <xsl:if test=\"position()&gt;1\">,</xsl:if>\n" + 
			"   <xsl:apply-templates select=\".\"/>\n" + 
			" </xsl:for-each>\n" + 
			"<xsl:text>\n" + 
			"</xsl:text>\n" + 
			" </xsl:if>\n" + 
			" \n" + 
			" <xsl:if test=\"intergenic_consequences\">\n" + 
			" <xsl:text>VEPREST_INTERGENIC_CSQ=</xsl:text>\n" + 
			" <xsl:for-each select=\"intergenic_consequences\">\n" + 
			"   <xsl:if test=\"position()&gt;1\">,</xsl:if>\n" + 
			"   <xsl:apply-templates select=\".\"/>\n" + 
			" </xsl:for-each>\n" + 
			"<xsl:text>\n" + 
			"</xsl:text>\n" + 
			" </xsl:if>\n" + 
			" <xsl:if test=\"regulatory_feature_consequences\">\n" + 
			" <xsl:text>VEPREST_REGULATORY_CSQ=</xsl:text>\n" + 
			" <xsl:for-each select=\"regulatory_feature_consequences\">\n" + 
			"   <xsl:if test=\"position()&gt;1\">,</xsl:if>\n" + 
			"   <xsl:apply-templates select=\".\"/>\n" + 
			" </xsl:for-each>\n" + 
			"<xsl:text>\n" + 
			"</xsl:text>\n" + 
			" </xsl:if>\n" + 
			" \n" + 
			"</xsl:template>\n" + 
			"\n" + 
			"<xsl:template match=\"transcript_consequences|intergenic_consequences|regulatory_feature_consequences\">\n" + 
			"<xsl:for-each select=\"@*\">\n" + 
			"   <xsl:if test=\"position()&gt;1\">|</xsl:if>\n" + 
			"   <xsl:value-of select=\"name(.)\"/>\n" + 
			"   <xsl:text>|</xsl:text>\n" + 
			"   <xsl:value-of select=\".\"/>\n" + 
			" </xsl:for-each>\n" + 
			" <xsl:if test=\"consequence_terms\">\n" + 
			" <xsl:text>|consequences|</xsl:text>\n" + 
			" <xsl:for-each select=\"consequence_terms\">\n" + 
			"   <xsl:if test=\"position()&gt;1\">&amp;</xsl:if>\n" + 
			"   <xsl:apply-templates select=\".\"/>\n" + 
			" </xsl:for-each>\n" + 
			" </xsl:if>\n" + 
			"</xsl:template>\n" + 
			"\n" + 
			"<xsl:template match=\"colocated_variants\">\n" + 
			"<xsl:for-each select=\"@*\">\n" + 
			"   <xsl:if test=\"position()&gt;1\">|</xsl:if>\n" + 
			"   <xsl:value-of select=\"name(.)\"/>\n" + 
			"   <xsl:text>|</xsl:text>\n" + 
			"   <xsl:value-of select=\".\"/>\n" + 
			" </xsl:for-each>\n" + 
			"</xsl:template>\n" + 
			"\n" + 
			"</xsl:stylesheet>\n"
			;
	private static final String DEFAULT_TAGS="VEPREST_TRANSCRIPT,VEPREST_COLOCATED_VARIANTS,VEPREST_INTERGENIC_CSQ,VEPREST_REGULATORY_CSQ";
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-s","--server"},description="Ensemble REST server url.")
	private String server = "http://grch37.rest.ensembl.org";

	@Parameter(names={"-e","--extension"},description="REST server URL Path extension")
	private String extension = "/vep/homo_sapiens/region";

	@Parameter(names={"-n","--batchSize"},description="batch size. How many variant to send in one HTTP query")
	private int batchSize = 100 ;

	@Parameter(names={"-x","--xslt","--template"},description="[20180214] XSLT stylesheet template that will be used to transform the XML response from ensembl to a set of lines."
			+ "Each line in the generated document must starts with a tag (TAG=...) defined with the --tag section. "
			+ "Empty lines or lines starting with a hash will be ignored. "
			+ "The XSL output method must be text. "
			+ "The processed node is a <data> node in the Rest response. "
			+ "Example at: https://gist.github.com/lindenb/12dea8e701d18280e8cb0b65270064f6 ."
			+ "The default XSL template is embedded in the code: "+
			"https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/tools/ensembl/VcfEnsemblVepRest.java"
			)
	private File xslTemplateFile = null;
	
	@Parameter(names={"-t","--tags","--tag"},description="[20180214] comma/space separated of tag that will be generated by the XSL stylesheet.")
	private  String tagListString  = DEFAULT_TAGS;
	
	@Parameter(names={"-T","--tee"},description="'Tee' xml response to stderr")
	private boolean teeResponse = false;
	@Parameter(names={"-nofail"},description="[20180213] Do not fail on network error")
	private boolean ignoreNetworkErrors = false;


	private DocumentBuilder documentBuilder;
	private Transformer xsltTransformer = null;
	private CloseableHttpClient httpClient = null;
	private final Set<String> outputTags = new HashSet<>();
	
	
	private class EnsVepPrediction
		{
		final String input;
		final Map<String,Set<String>> tag2infoLines;
		EnsVepPrediction(final String input)
			{
			this.input = input;
			this.tag2infoLines = new HashMap<>(VcfEnsemblVepRest.this.outputTags.size());
			for(final String tag:VcfEnsemblVepRest.this.outputTags)
				{
				this.tag2infoLines.put(tag, new HashSet<>());
				}
			}	
		}
	
	
	
	private Map<String,EnsVepPrediction> parseVepPredictions(final Document dom) {
		final Map<String,EnsVepPrediction> preds = new HashMap<>();
		final Element root= dom.getDocumentElement();
		if(root==null) {
			if(ignoreNetworkErrors) {
				LOG.error("empty dom");
				return Collections.emptyMap();
				}
			throw new IllegalStateException("empty dom");
			}
		if(!root.getNodeName().equals("opt")) {
			final String msg = "root is not <opt> but <"+root.getNodeName()+">";
			if(ignoreNetworkErrors) {
				LOG.error( msg);
				return Collections.emptyMap();
				}
			throw new IllegalStateException( msg);
		}
		for(Node c1 = root.getFirstChild();c1!=null;c1=c1.getNextSibling())
			{
			if(c1.getNodeType()!=Node.ELEMENT_NODE) continue;
			final Element e1 = Element.class.cast(c1);
			if(e1.getNodeName().equals("data"))
				{
				if(!e1.hasAttribute("input"))
					{
					LOG.warning("found <data> without @input");
					continue;
					}	
				try
					{
					final EnsVepPrediction pred = new EnsVepPrediction(e1.getAttribute("input"));

					final StringWriter sw = new StringWriter();
					final StreamResult result = new StreamResult(sw);
					this.xsltTransformer.transform(new DOMSource(e1), result);
					sw.close();
					final BufferedReader br= new BufferedReader(new StringReader(sw.toString()));
					String line;
					while((line=br.readLine())!=null) {
						if(StringUtil.isBlank(line) || line.startsWith("#")) continue;
						final int eq = line.indexOf("=");
						if(eq==-1 ) throw new TransformerException("Cannot find '=' in "+line);
						final String tag = line.substring(0,eq);
						if(!pred.tag2infoLines.containsKey(tag))
							{
							throw new TransformerException("unknown tag '"+tag+"' in "+line+". Defined are "+pred.tag2infoLines.keySet());
							}
						final String value=line.substring(eq+1).trim();
						if(StringUtil.isBlank(value)) continue;
						if(value.contains(" ") || value.contains("\t")) {
							throw new TransformerException("spaces in value: '"+line+"'.");
							}
						pred.tag2infoLines.get(tag).add(value);
						}
					br.close();
					if(preds.containsKey(pred.input))
						{
						LOG.error("duplicate data/@input: "+pred.input+" in XML");
						}
					preds.put(pred.input,pred);
					}
				catch(final Exception err)
					{
					final String msg= "XSLT Transformation failed : "+err.getMessage();
					if(ignoreNetworkErrors) {
						LOG.error( msg,err);
						continue;
						}
					throw new IllegalStateException(msg,err);
					}
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
			append(ctx.getReference().getDisplayString()).
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
				sb.append(alts.get(j).getDisplayString());
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
	
	private Map<String,EnsVepPrediction> vep(final List<VariantContext> contexts) throws IOException
		{
		return parseVepPredictions(callVepToDom(contexts,false));
		}
	
	@Override
	public int doWork(final List<String> args) {
	try {
		this.outputTags.addAll(Arrays.stream(this.tagListString.split("[, \t;]+")).filter(S->!StringUtil.isBlank(S)).collect(Collectors.toSet()));
		for(final String tag:this.outputTags)
			{
			if(!tag.matches("[A-Za-z][A-Za-z0-9_]*"))
				{
				LOG.error("Bad tag name "+tag);
				return -1;
				}
			}
		if(this.outputTags.isEmpty()) {
			LOG.error("no output tag defined");
			return -1;
		}
			
		
		final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		this.documentBuilder=dbf.newDocumentBuilder();		
		final  TransformerFactory trf=TransformerFactory.newInstance();
		final Source stylesheetSource;
		if(this.xslTemplateFile!=null)
			{
			stylesheetSource = new StreamSource(this.xslTemplateFile); 
			}
		else
			{
			stylesheetSource = new StreamSource(new StringReader(DEFAULT_XSLT_TEMPLATE)); 
			}
		final Templates templates = trf.newTemplates(stylesheetSource);
		this.xsltTransformer = templates.newTransformer();
		
		/** create http client */
		this.httpClient = HttpClients.createSystem();//createDefault();
		
		return doVcfToVcf(args, this.outputFile);
		}
	catch(final Exception err) {
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
	protected int doVcfToVcf(final String inputName, final VCFIterator vcfIn, final VariantContextWriter out) {
	    try {
		final VCFHeader header=vcfIn.getHeader();
		final List<VariantContext> buffer=new ArrayList<>(this.batchSize+1);
		final VCFHeader h2= new VCFHeader(header);
		addMetaData(h2);
		for(final String tag:this.outputTags)
			{
			h2.addMetaDataLine(new VCFInfoHeaderLine(
					tag,
					VCFHeaderLineCount.UNBOUNDED,
					VCFHeaderLineType.String,
					"__CUSTOM_DESCRIPTION__"+tag
					));
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
					final Map<String,EnsVepPrediction> input2pred = vep(buffer);
					for(final VariantContext ctx2:buffer)
						{
						final String inputStr = createInputContext(ctx2);
						final EnsVepPrediction pred= input2pred.get(inputStr);
						
						final VariantContextBuilder vcb=new VariantContextBuilder(ctx2);
						for(final String tag:this.outputTags)
							{
							vcb.rmAttribute(tag);
							}
						
						if(pred==null)
							{
							LOG.info("No Annotation found for "+inputStr);
							out.add(vcb.make());
							continue;
							}
						
						for(final String tag:this.outputTags)
							{
							final Set<String> info = pred.tag2infoLines.get(tag);
							if(info==null || info.isEmpty()) continue;
							vcb.attribute(tag, new ArrayList<>(info));
							}
						
						out.add(vcb.make());
						} // end of loop over variants
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
	
	
	public static void main(final String[] args) {
		new VcfEnsemblVepRest().instanceMainWithExit(args);
	}
	}
