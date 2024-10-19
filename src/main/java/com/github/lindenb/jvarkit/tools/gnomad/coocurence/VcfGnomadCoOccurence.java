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
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.concurrent.TimeUnit;

import org.apache.http.HttpResponse;
import org.apache.http.client.HttpClient;
import org.apache.http.client.methods.HttpPost;
import org.apache.http.entity.StringEntity;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.MathUtils;
import com.github.lindenb.jvarkit.net.ContentType;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.internal.Streams;
import com.google.gson.stream.JsonWriter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StopWatch;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC
 

## Example

```
java -jar dist/jvarkit.jarvcfgnomadcoocurence \
 src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz \
 --bed in.bed

(...)
    {
        "data": {
            "variant_cooccurrence": {
                "variant_ids": [
                    "1-905958-G-T",
                    "1-905962-C-A"
                ],
                "genotype_counts": [
                    124777,
                    1,
                    0,
                    8,
                    0,
                    0,
                    0,
                    0,
                    0
                ],
                "haplotype_counts": [
                    249563,
                    1,
                    8,
                    0
                ],
                "p_compound_heterozygous": 1,
                "populations": [
                    {
                        "id": "afr",
                        "genotype_counts": [
                            8042,
                            0,
                            0,
                            8,
                            0,
                            0,
                            0,
                            0,
                            0
                        ],
                        "haplotype_counts": [
                            16092,
                            0,
                            8,
                            0
                        ],
                        "p_compound_heterozygous": null
                    },
(...)

```

END_DOC
 */
@Program(name="vcfgnomadcoocurence",
	description="Programmatic use of gnomad co-ocurence",
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
	@Parameter(names={"-B","--bed","--regions"},description="Bed intervals to scan",required=true)
	protected Path bedPath=null;
	@Parameter(names={"--sleep"},description="Sleep 'x' seconds between each call")
	protected int seconds =1;
	@Parameter(names={"--sleep2"},description="Sleep 'x' seconds on error")
	protected int seconds2 =120;

	
	private static String toString(final VariantContext ctx) {
		String ctg=ctx.getContig();
		if(ctg.startsWith("chr")) ctg=ctg.substring(3);
		return String.join("-", ctg,
			String.valueOf(ctx.getStart()),
			ctx.getAlleles().get(0).getBaseString(),
			ctx.getAlleles().get(1).getBaseString()
			);
		}
	
	private OptionalDouble report_p_compound_heterozygous(final JsonElement root) {
		if(root.isJsonArray()) {
			final JsonArray a = root.getAsJsonArray();
			for(int i=0;i< a.size();++i) {
				OptionalDouble d=report_p_compound_heterozygous(a.get(i));
				if(d.isPresent()) return d;
			}
		}
		else if(root.isJsonArray()) {
			final JsonObject o = root.getAsJsonObject();
			for(Map.Entry<String, JsonElement> kv:o.entrySet()) {
				if(kv.getKey().equals("p_compound_heterozygous")) {
					JsonElement v=kv.getValue();
					if(!v.isJsonNull() && v.isJsonPrimitive()) {
						Number f = v.getAsNumber();
						if(f.doubleValue()!=1.0) {
							return OptionalDouble.of(f.doubleValue());
							}
						}
					}
				OptionalDouble d= report_p_compound_heterozygous(kv.getValue());
				if(d.isPresent()) return d;
			}
		}
	return OptionalDouble.empty();
	}
	
	private boolean fetch(final HttpClient client,final JsonWriter jsw,final VariantContext ctx1,final VariantContext ctx2) {
		final int try_count=10;
		for(int n_try=0;n_try<try_count;++n_try) {
			if(n_try>0)  {
				LOG.error("ntry "+(n_try+1)+"/"+try_count);
				sleep(seconds2);
				}
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
					/*
					.append("    exome {\n")
					.append("      ac\n")
					.append("      an\n")
					.append("    }\n") */
					.append("    genome {\n")
					.append("      ac\n")
					.append("      an\n")
					.append("    }\n")
					/*
					.append("    multi_nucleotide_variants {\n")
					.append("      combined_variant_id\n")
					.append("      other_constituent_snvs\n")
					.append("    }\n")
					.append("    transcript_consequences {\n")
					.append("      gene_id\n")
					.append("      gene_version\n")
					.append("      gene_symbol\n")
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
					.append("      transcript_version\n")
					.append("    }\n")
					*/
					.toString();


			final StringBuilder w=new StringBuilder();
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


				final StringEntity params = new StringEntity(query.toString());
				request.setEntity(params);

				final HttpResponse response = client.execute(request);
				if(response.getStatusLine().getStatusCode()!=200) {
					LOG.error("BAD RESPONSE:\n\t"+response.getStatusLine());
					
					continue;
				}
				try (InputStream in = response.getEntity().getContent()) {
					try (InputStreamReader ir = new InputStreamReader(in)) {
						final JsonParser parser = new JsonParser();

						final JsonElement root;

						try {
							root = parser.parse(ir);
						} catch (Throwable err0) {
							LOG.error(err0);
							err0.printStackTrace();
							continue;
						}
						if (!root.isJsonObject() || !root.getAsJsonObject().has("data")) {
							LOG.warning("No data found for " + vc1 + " " + vc2 + "\n\t" + root);
						} else {
							final OptionalDouble p_compound_heterozygous = report_p_compound_heterozygous(root);
							if(p_compound_heterozygous.isPresent()) {
								LOG.info("Interesting ! "+vc1+" "+vc2+" "+p_compound_heterozygous);
								}
							
							
							final JsonObject data = root.getAsJsonObject().get("data").getAsJsonObject();
							for (int i = 1; i <= 2; i++) {
								if (!data.has("variant" + i)) {
									continue;
									}
								final JsonObject o = data.get("variant" + i).getAsJsonObject();
								final VariantContext ctx = (i == 1 ? ctx1 : ctx2);
								o.addProperty("location",toString(ctx));
								o.addProperty("contig", ctx.getContig());
								o.addProperty("start", ctx.getStart());
								o.addProperty("end", ctx.getEnd());
								o.addProperty("ref", ctx.getAlleles().get(0).getBaseString());
								o.addProperty("alt", ctx.getAlleles().get(1).getBaseString());
							}
						}

						Streams.write(root, jsw);
						jsw.flush();
					}
				}
				return true;
			} catch(Throwable err) {
				LOG.error(err);
			}
		} 
	return false;
	}
	
	private static void sleep(int seconds) {
		if(seconds>0) {
			try {
				TimeUnit.SECONDS.sleep(seconds);
				}
			catch(Throwable err) {
				
			}
		}
	}
	
	
	
	@Override
	public int doWork(final List<String> args) {
			
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
				final List<Locatable> intervals = new ArrayList<>();
				try(BedLineReader blr = new BedLineReader(this.bedPath)) {
					blr.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
					blr.forEachRemaining(B->intervals.add(B));
					}
				
				if(intervals.isEmpty()) {
					LOG.warn("no interval in "+this.bedPath);
					}
				
				
				int status = 0;
				/** create http client */
				try(CloseableHttpClient httpClient = HttpClients.createSystem()) {
					try(PrintWriter pw =super.openPathOrStdoutAsPrintWriter(outputFile)) {
						try(JsonWriter jsw = new JsonWriter(pw)) {
							jsw.beginArray();
							
							for(final Locatable interval:intervals) {								
								final List<VariantContext> variants = new ArrayList<>();
							
								try(CloseableIterator<VariantContext> iter=reader.query(interval)) {
									while(iter.hasNext()) {
										final VariantContext ctx = iter.next();
										if(!AcidNucleics.isATGC(ctx.getReference())) continue;
										if(ctx.getReference().length()!=1) continue;
										for(Allele alt:ctx.getAlternateAlleles()) {
											if(!AcidNucleics.isATGC(alt)) continue;
											if(alt.length()!=1) continue;

											final VariantContextBuilder vcb = new VariantContextBuilder(
												input,
												ctx.getContig(),
												ctx.getEnd(),
												ctx.getEnd(),
												Arrays.asList(ctx.getReference(),alt));
											variants.add(vcb.make());
											}
										}
									}
								final StopWatch stopWatch = new StopWatch();
								stopWatch.start();
								final long n_expect=variants.size()<2?0L:MathUtils.combination(variants.size(),2).longValue();
								long n_count=0L;
								for(int x=0;x+1<variants.size() && status==0;++x ) {
									final VariantContext ctx1 = variants.get(x);
									for(int y=x+1;y<variants.size() && status==0;++y ) {
										final VariantContext ctx2 = variants.get(y);
										if(ctx1.getStart()==ctx2.getStart()) {
											++n_count;
											continue;
											}
										LOG.info(toString(ctx1)+"/"+toString(ctx2)+" "+StringUtils.niceInt(n_count)+"/"+StringUtils.niceInt(n_expect)+
												" remains "+(n_count==0?"xxxxx":""+(int)(((n_expect-n_count)*(stopWatch.getElapsedTimeSecs()/(double)n_count))))+" secs");
										if(!fetch(httpClient,jsw,ctx1,ctx2)) {
											LOG.error("Error");
											status=-1;
											break;
											}
										sleep(seconds);
										++n_count;
										}
									}
								stopWatch.stop();
								if(status!=0) break;
								} /* end loop over intervals */
							
							jsw.endArray();
							jsw.flush();
							}
						pw.flush();
						}
					}
				return status;
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
