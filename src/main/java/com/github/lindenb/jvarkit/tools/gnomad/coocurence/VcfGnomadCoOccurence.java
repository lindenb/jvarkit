/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.gnomad.coocurence;

import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.OptionalDouble;
import java.util.function.Function;

import org.apache.http.HttpResponse;
import org.apache.http.client.HttpClient;
import org.apache.http.client.methods.HttpPost;
import org.apache.http.entity.StringEntity;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.net.ContentType;
import com.github.lindenb.jvarkit.samtools.util.LocatableUtils;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC
 

`

END_DOC
 */
@Program(name="vcfgnomadcoocurence",
	description="Programmatic use of gnomad coocurence",
	keywords={"vcf","gnomad"},
	modificationDate="20241018",
	creationDate="20241018",
	jvarkit_amalgamion =  true,
	menu="VCF Manipulation"
)
public class VcfGnomadCoOccurence extends Launcher {
	private static final Logger LOG = Logger.build(VcfGnomadCoOccurence.class).make();

	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	protected Path outputFile=null;
	@Parameter(names={"-r","--regions","--region"},description="Intervals: one interval",required=true)
	protected String intervalStr="";
	
	
	static class VariantCoOccurence {
		List<Integer> genotype_counts;
		List<Integer> haplotype_counts;
		OptionalDouble p_compound_heterozygous=OptionalDouble.empty();
	}
	
	private static String toString(VariantContext ctx) {
		String ctg=ctx.getContig();
		if(ctg.startsWith("chr")) ctg=ctg.substring(3);
		return String.join("-", ctg,
			String.valueOf(ctx.getStart()),
			ctx.getAlleles().get(0).getBaseString(),
			ctx.getAlleles().get(1).getBaseString()
			);
		}
	
	private static <T> List<T> toList(JsonArray array,Function<JsonElement,T> converter ) {
		List<T> L= new ArrayList<>(array.size());
		for(int i=0;i< array.size();i++) {
			L.add(converter.apply(array.get(i)));
			}
		return L;
		}
	
	private void fetch(final HttpClient client,VariantContext ctx1,VariantContext ctx2) {
		final String vc1 = toString(ctx1);
		final String vc2 = toString(ctx2);
		final HttpPost request = new HttpPost("https://gnomad.broadinstitute.org/api/");
		request.addHeader("User-Agent", IOUtils.getDefaultUserAgent());
		request.addHeader("Accept"," */*");
		request.addHeader("Accept-Language"," en-US,en;q=0.5");
		request.addHeader("Accept-Encoding"," gzip, deflate, br, zstd");
		request.addHeader("Referer","https://gnomad.broadinstitute.org/variant-cooccurrence?dataset=gnomad_r2_1&variant="+
				vc1+
				"&variant="+
				vc2
				);
		request.addHeader("Content-Type",ContentType.APPLICATION_JSON.getMimeType());
		request.addHeader("Origin","https://gnomad.broadinstitute.org");
		request.addHeader("DNT","1");
		request.addHeader("Sec-GPC","1");
		request.addHeader("Connection"," keep-alive");
		request.addHeader("Sec-Fetch-Dest","empty");
		request.addHeader("Sec-Fetch-Mode","cors");
		request.addHeader("Sec-Fetch-Site","same-origin");
		request.addHeader("Priority","u=0");
		request.addHeader("TE","trailers");

		final String components = new StringBuilder()
	/*	.append("    exome {\n")
		.append("      ac\n")
		.append("      an\n")
		.append("    }\n")
		.append("    genome {\n")
		.append("      ac\n")
		.append("      an\n")
		.append("    }\n")
		.append("    multi_nucleotide_variants {\n")
		.append("      combined_variant_id\n")
		.append("      other_constituent_snvs\n")
		.append("    }\n")*/
		.append("    transcript_consequences {\n")
		.append("      gene_id\n")
		.append("      gene_version\n")
		/*.append("      gene_symbol\n")
		.append("      hgvs\n")
		.append("      hgvsc\n")
		.append("      hgvsp\n")
		.append("      is_canonical\n")
		.append("      is_mane_select\n")
		.append("      is_mane_select_version\n")
		.append("      lof\n")
		.append("      lof_flags\n")
		.append("      lof_filter\n")
		.append("      major_consequence\n")
		.append("      polyphen_prediction\n")
		.append("      sift_prediction\n")
		.append("      transcript_id\n")
		.append("      transcript_version\n")*/
		.append("    }\n")
		.toString();
		
		
		StringBuilder w=new StringBuilder();
		
		w.append("query VariantCooccurrence($variants: [String!]!, $variant1: String!, $variant2: String, $datasetId: DatasetId!) {\n");
		w.append("  variant_cooccurrence(variants: $variants, dataset: $datasetId) {\n");
		w.append("    variant_ids\n");
		w.append("    genotype_counts\n");
		w.append("    haplotype_counts\n");
		w.append("    p_compound_heterozygous\n");
		w.append("    populations {\n");
		w.append("      id\n");
		w.append("      genotype_counts\n");
		w.append("      haplotype_counts\n");
		w.append("      p_compound_heterozygous\n");
		w.append("    }\n");
		w.append("  }\n");
		w.append("  variant1: variant(variantId: $variant1, dataset: $datasetId) {\n");
		w.append(components).append("\n");
		w.append("  }\n");
		w.append("  variant2: variant(variantId: $variant2, dataset: $datasetId) {\n");
		w.append(components);
		w.append("  }\n");
		w.append("}\n");

		//String q2="\nquery VariantCooccurrence($variants: [String!]!, $variant1: String!, $variant2: String, $datasetId: DatasetId!) {\n  variant_cooccurrence(variants: $variants, dataset: $datasetId) {\n    variant_ids\n    genotype_counts\n    haplotype_counts\n    p_compound_heterozygous\n    populations {\n      id\n      genotype_counts\n      haplotype_counts\n      p_compound_heterozygous\n    }\n  }\n  variant1: variant(variantId: $variant1, dataset: $datasetId) {\n    exome {\n      ac\n      an\n    }\n    genome {\n      ac\n      an\n    }\n    multi_nucleotide_variants {\n      combined_variant_id\n      other_constituent_snvs\n    }\n    transcript_consequences {\n      gene_id\n      gene_version\n      gene_symbol\n      hgvs\n      hgvsc\n      hgvsp\n      is_canonical\n      is_mane_select\n      is_mane_select_version\n      lof\n      lof_flags\n      lof_filter\n      major_consequence\n      polyphen_prediction\n      sift_prediction\n      transcript_id\n      transcript_version\n    }\n  }\n  variant2: variant(variantId: $variant2, dataset: $datasetId) {\n    exome {\n      ac\n      an\n    }\n    genome {\n      ac\n      an\n    }\n    multi_nucleotide_variants {\n      combined_variant_id\n      other_constituent_snvs\n    }\n    transcript_consequences {\n      gene_id\n      gene_version\n      gene_symbol\n      hgvs\n      hgvsc\n      hgvsp\n      is_canonical\n      is_mane_select\n      is_mane_select_version\n      lof\n      lof_flags\n      lof_filter\n      major_consequence\n      polyphen_prediction\n      sift_prediction\n      transcript_id\n      transcript_version\n    }\n  }\n}\n";
		
		JsonObject query=new JsonObject();
		query.addProperty("operationName", "VariantCooccurrence");
		query.addProperty("query",w.toString());
		
		final JsonArray variants = new JsonArray();
		variants.add(vc1);
		variants.add(vc2);
		
		final JsonObject variables = new JsonObject();
		variables.add("variants",variants);
		variables.addProperty("variant1",vc1);
		variables.addProperty("variant2",vc2);
		variables.addProperty("datasetId","gnomad_r2_1");
		query.add("variables",variables);
		
		try {
			
			JsonParser parser = new JsonParser();
			parser.parse(query.toString());
			
		 StringEntity params = new StringEntity(query.toString());
		 request.setEntity(params);
		 
		  HttpResponse response = client.execute(request);
		  try(InputStream in=response.getEntity().getContent()) {
			  try(InputStreamReader ir=new InputStreamReader(in)) {
			  JsonElement je=parser.parse(ir);
			  JsonObject obj1 = je.getAsJsonObject();
			  if(!obj1.has("data")) {
				  
			  	}
			  else
				  {
				  JsonObject data = obj1.get("data").getAsJsonObject();
				  List<Integer> genotype_counts  = toList(data.get("genotype_counts").getAsJsonArray(),O->O.getAsInt());
				  List<Integer> haplotype_counts  = toList(data.get("haplotype_counts").getAsJsonArray(),O->O.getAsInt());
				  }
			  }
		  }
		} catch(Throwable err) {
			LOG.error(err);
			err.printStackTrace();
		}
		 
		
	}
	
	@Override
	public int doWork(List<String> args) {
			final String input = super.oneAndOnlyOneFile(args);
			try(final VCFReader reader= VCFReaderFactory.makeDefault().open(Paths.get(input), true)) {
				final VCFHeader header = reader.getHeader();
				final SAMSequenceDictionary dict = header.getSequenceDictionary();
				if(dict!=null ) {
					if(!SequenceDictionaryUtils.isGRCh37(dict)) {
						LOG.info("doesn't look like a GRCH37 VCF file");
						return -1;
						}
					}

				final Locatable interval = LocatableUtils.parse(this.intervalStr,dict);
				
				final List<VariantContext> variants = new ArrayList<>();
			
				try(CloseableIterator<VariantContext> iter=reader.query(interval)) {
					while(iter.hasNext()) {
						VariantContext ctx = iter.next();
						if(!AcidNucleics.isATGC(ctx.getReference())) continue;
						if(ctx.getReference().length()!=1) continue;
						for(Allele alt:ctx.getAlternateAlleles()) {
							if(!AcidNucleics.isATGC(alt)) continue;
							if(alt.length()!=1) continue;

							VariantContextBuilder vcb = new VariantContextBuilder(
									input,
									ctx.getContig(),
									ctx.getEnd(),
									ctx.getEnd(),
									Arrays.asList(ctx.getReference(),alt));
							variants.add(vcb.make());
							}
						}
					}
				
				/** create http client */
				try(CloseableHttpClient httpClient = HttpClients.createSystem()) {
					for(int x=0;x+1<variants.size();++x ) {
						VariantContext ctx1 = variants.get(x);
						for(int y=x+1;y<variants.size();++y ) {
							VariantContext ctx2 = variants.get(y);
							if(ctx1.getStart()==ctx2.getStart()) continue;
							fetch(httpClient,ctx1,ctx2);
							}
						}
					}
				
				
				
				return 0;
				}
			catch(Throwable err) {
				LOG.error(err);
				return -1;
				}
			
		
		
		
		}
	
	
	
	
	
	

public static void main(final String[] args) {
	new VcfGnomadCoOccurence().instanceMainWithExit(args);
	}
}
