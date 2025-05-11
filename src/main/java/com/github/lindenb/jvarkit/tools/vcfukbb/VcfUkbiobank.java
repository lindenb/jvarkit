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
package com.github.lindenb.jvarkit.tools.vcfukbb;

import java.net.URI;
import java.net.http.HttpClient;
/**
BEGIN_DOC



END_DOC
 */
import java.net.http.HttpClient.Version;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
/**
BEGIN_DOC

# WARNING

I cannot explain why this software doesn't run with jdk17 (I always get a http error 403) but runs with jdk23...

The server might not like too many requests. Use a your own risk.

## Example:

```
$ cat jeter.vcf | grep -v "#"
chr22	10510356	.	T	A	.	.	.
chr22	10526438	.	A	G	.	.	.
```

```
java -jar jvarkit.jar vcfukbb jeter.vcf | grep -v "#"
chr22	10510356	.	T	A	.	.	UKBB_AFR_AC=288;UKBB_AFR_AF=0.319;UKBB_ALL_AC=12606;UKBB_ALL_AF=0.915;UKBB_ASJ_AC=71;UKBB_ASJ_AF=0.934;UKBB_EAS_AC=779;UKBB_EAS_AF=0.986;UKBB_NFE_AC=9989;UKBB_NFE_AF=0.954;UKBB_OTH_AC=448;UKBB_OTH_AF=0.907;UKBB_SAS_AC=1031;UKBB_SAS_AF=0.982
chr22	10526438	.	A	G	.	.	UKBB_AFR_AC=13052;UKBB_AFR_AF=0.863;UKBB_ALL_AC=536812;UKBB_ALL_AF=0.986;UKBB_ASJ_AC=3261;UKBB_ASJ_AF=0.964;UKBB_EAS_AC=3119;UKBB_EAS_AF=0.960;UKBB_NFE_AC=495657;UKBB_NFE_AF=0.990;UKBB_OTH_AC=9134;UKBB_OTH_AF=0.964;UKBB_SAS_AC=12589;UKBB_SAS_AF=0.991
```

##

END_DOC

 */
@Program(
		name="vcfukbb",
		description="annotates an VCF with the https://afb.ukbiobank.ac.uk/ ukbiobank server. The server might not like too many requests. Use a your own risk. Doesn't work with jdk17 (?!)",
		keywords={"ukbiobank","vcf"},
		creationDate="20250424",
		modificationDate="20250424",
		generate_doc = true,
		jvarkit_amalgamion = true
		)
public class VcfUkbiobank extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.of(VcfUkbiobank.class);
	private final HttpClient httpClient = HttpClient.newHttpClient();
	private final JsonParser jsonparser = new JsonParser();

	@Parameter(names = { "--skip-filtered" }, description = "Skip FILTERed variants")
	private boolean skip_filtered_variants = false;
	@Parameter(names = { "--api" }, description = "API base url")
	private String api_url = "https://afb.ukbiobank.ac.uk/api";
	@Parameter(names = { "--seconds" }, description = "wait 'x' seconds between each call to the API")
	private int max_seconds = 10;
	@Parameter(names = { "--retry" }, description = "max-retry")
	private int max_retries = 10;
	@Parameter(names = { "--debug" }, description = "debug")
	private boolean debug = false;

	// private Locatable bufferInterval = null;
	// private final List<VariantContextBuilder> buffer=new ArrayList<>();

	private static class PopulationData {
		final String popName;
		final VCFInfoHeaderLine allele_count_info;
		int[] allele_count;
		final VCFInfoHeaderLine allele_num_info;
		Integer allele_num;
		final VCFInfoHeaderLine allele_freq_info;
		double[] allele_freq;
		final VCFInfoHeaderLine nHomozygotes_info;
		Integer nHomozygotes;
		final VCFInfoHeaderLine nHemiAlt_info;
		Integer nHemiAlt;

		PopulationData(String popName, String tag) {
			this.popName = popName;
			this.allele_count_info = new VCFInfoHeaderLine("UKBB_" + tag + "_AC", VCFHeaderLineCount.A,
					VCFHeaderLineType.Integer, "Allele Count for " + popName);
			this.allele_num_info = new VCFInfoHeaderLine("UKBB_" + tag + "_AN", 1, VCFHeaderLineType.Integer,
					"Allele Number for " + popName);
			this.allele_freq_info = new VCFInfoHeaderLine("UKBB_" + tag + "_AF", VCFHeaderLineCount.A,
					VCFHeaderLineType.Integer, "Allele Frequency for " + popName);
			this.nHomozygotes_info = new VCFInfoHeaderLine("UKBB_" + tag + "_HOM", 1, VCFHeaderLineType.Integer,
					"Number of homozygotes for " + popName);
			this.nHemiAlt_info = new VCFInfoHeaderLine("UKBB_" + tag + "_HEMI", 1, VCFHeaderLineType.Integer,
					"Number of hemizygote ALT for " + popName);
		}

		void reset(int n_alts) {
			allele_count = new int[n_alts];
			Arrays.fill(allele_count, 0);
			allele_freq = new double[n_alts];
			Arrays.fill(allele_freq, 0);
			allele_num = null;
			nHomozygotes = null;
		}

	}

	private static class Response {
		// int responseCode;
		// String errorMsg = "";
		JsonElement root = null;
		boolean fail_flag = false;
	}

	private void build(PopulationData pop, VariantContextBuilder vcb) {
		vcb.attribute(pop.allele_count_info.getID(), pop.allele_count);
		vcb.attribute(pop.allele_freq_info.getID(), pop.allele_freq);
		if (pop.allele_num != null)
			vcb.attribute(pop.allele_num_info.getID(), pop.allele_num);
	}

	private void populationDetails(JsonObject root, int alt_idx, PopulationData pop) {
		if (root.has("alleleCount") && root.get("alleleCount").isJsonPrimitive()) {
			pop.allele_count[alt_idx] = root.get("alleleCount").getAsInt();
		}

		if (root.has("alleleFreq") && root.get("alleleFreq").isJsonPrimitive()) {
			pop.allele_freq[alt_idx] = root.get("alleleFreq").getAsDouble();
		}

		if (pop.allele_num != null && root.has("alleleNum") && root.get("alleleNum").isJsonPrimitive()) {
			pop.allele_num = root.get("alleleNum").getAsInt();
		}
		if (pop.nHomozygotes != null && root.has("nHomozygotes")
				&& root.get("nHomozygotes").isJsonPrimitive()) {
			pop.nHomozygotes = root.get("nHomozygotes").getAsInt();
		}
		if (pop.nHemiAlt != null && root.has("nHemiAlt") && root.get("nHemiAlt").isJsonPrimitive()) {
			pop.nHemiAlt = root.get("nHemiAlt").getAsInt();
		}

	}

	private void populationDetails(final JsonElement root, int alt_idx,
			final Map<String, PopulationData> hash) {
		if (root == null || !root.isJsonArray())
			return;
		final JsonArray array = root.getAsJsonArray();
		for (int i = 0; i < array.size(); ++i) {
			if (!array.get(i).isJsonObject())
				continue;
			com.google.gson.JsonObject o2 = array.get(i).getAsJsonObject();
			if (!o2.has("population") || !o2.get("population").isJsonPrimitive())
				continue;
			final String pop = o2.get("population").getAsString();
			if (!hash.containsKey(pop))
				continue;
			populationDetails(o2, alt_idx, hash.get(pop));
		}
	}

	private String getContig(final VariantContext ctx) {
		return ctx.getContig().startsWith("chr") ? ctx.getContig() : "chr" + ctx.getContig();
	}

	/*
	 * private boolean hasVariant(final VariantContext ctx,final Allele alt) {
	 * if(bufferInterval==null || !bufferInterval.overlaps(ctx)) { bufferInterval =
	 * new SimpleInterval(ctx.get) } }
	 */

	private Response query(final VariantContext ctx, final Allele alt) {
		final Response apiReturn = new Response();
		final String contig = getContig(ctx);

		final String query = "{\"query\":\"\\nquery VariantSplitPopulationQuery($chrom: String!, $pos: Int!, $ref: String!, $alt: String!) {\\n    variant(chrom: $chrom, pos: $pos, ref: $ref, alt: $alt) {\\n      Chrom\\n      Pos\\n      Ref\\n      Alt\\n      alleleCount\\n      alleleNum\\n      alleleFreq\\n      nHomozygotes\\n      maxImpact\\n      maxConsequence\\n      HGVSp\\n      geneSymbol\\n      populationDetails {\\n        population\\n        alleleCount\\n        alleleNum\\n        alleleFreq\\n        nHomozygotes\\n        nHemiAlt\\n      }\\n    }\\n  } \\n\",\"variables\":{\"chrom\":\""
				+ contig + "\",\"pos\":" + ctx.getStart() + ",\"ref\":\"" + ctx.getReference().getBaseString()
				+ "\",\"alt\":\"" + alt.getBaseString() + "\"}}";
		if (debug)
			LOG.info(query);
		this.jsonparser.parse(query);// check ok

		{

			try {
				final HttpRequest.Builder hb = HttpRequest.newBuilder();
				hb.version(Version.HTTP_2);
				hb.header("User-Agent", IOUtils.getDefaultUserAgent());
				hb.header("Accept", "application/json, text/plain, */*");
				hb.header("Accept-Language", "en-US,en;q=0.5");
				hb.header("Content-Type", "application/json");
				hb.header("Origin", "https://afb.ukbiobank.ac.uk");
				hb.header("Referer", "https://afb.ukbiobank.ac.uk/variant/" + contig + "-" + ctx.getStart()
						+ "-" + ctx.getReference().getBaseString() + "-" + alt.getBaseString());
				hb.uri(URI.create(this.api_url));
				hb.POST(HttpRequest.BodyPublishers.ofString(query, StandardCharsets.UTF_8));
				final HttpRequest request = hb.build();
				final HttpResponse<String> resp = httpClient.send(request,
						HttpResponse.BodyHandlers.ofString());

				if (resp.statusCode() == 200) {
					final String jsonStr = resp.body();
					if (jsonStr.startsWith("<")) {
						apiReturn.fail_flag = true;
					}
					apiReturn.root = this.jsonparser.parse(jsonStr);
				} else {
					switch (resp.statusCode()) {
					case 403 /* forbidden */:
						apiReturn.fail_flag = true;
						break;
					default:
						apiReturn.fail_flag = true;
						break;
					}
					LOG.info("http returns " + resp.statusCode());
				}

			} catch (Throwable e) {
				if (debug)
					LOG.error(e);
			}
		}

		return apiReturn;
	}

	@Override
	protected int doVcfToVcf(final String inputName, final VCFIterator iter, final VariantContextWriter out) {
		final String GENERAL_POP = "All";
		final Map<String, PopulationData> popHash = new HashMap<>();
		for (PopulationData pop : Arrays.asList(new PopulationData(GENERAL_POP, "ALL"),
				new PopulationData("African", "AFR"), new PopulationData("Ashkenazi Jewish", "ASJ"),
				new PopulationData("East Asian", "EAS"), new PopulationData("Non-Finnish European", "NFE"),
				new PopulationData("Other", "OTH"), new PopulationData("South Asian", "SAS"))) {
			popHash.put(pop.popName, pop);
		}

		final VCFHeader header = iter.getHeader();
		long last_call = 0L;
		final SAMSequenceDictionary dict = header.getSequenceDictionary();
		if (dict != null && !SequenceDictionaryUtils.isGRCh38(dict)) {
			LOG.error("sequence dictionary doesn't look like hg38 ?!!:"
					+ SequenceDictionaryUtils.getBuildName(dict).orElse("undefined"));
			return -1;
		}
		final VCFHeader h2 = new VCFHeader(header);
		final VCFInfoHeaderLine info_ukbbstatus = new VCFInfoHeaderLine("UKBB_FAILED", 0,
				VCFHeaderLineType.Flag, "API Call Failed");
		h2.addMetaDataLine(info_ukbbstatus);

		for (PopulationData pop : popHash.values()) {
			h2.addMetaDataLine(pop.allele_count_info);
			h2.addMetaDataLine(pop.allele_num_info);
			h2.addMetaDataLine(pop.allele_freq_info);
			h2.addMetaDataLine(pop.nHomozygotes_info);
			h2.addMetaDataLine(pop.nHemiAlt_info);
		}

		JVarkitVersion.getInstance().addMetaData(this, h2);
		out.writeHeader(h2);
		while (iter.hasNext()) {
			final VariantContext ctx = iter.next();
			if (this.skip_filtered_variants && ctx.isFiltered()) {
				out.add(ctx);
				continue;
			}

			if (!ctx.getContig().matches("(chr)?[XY0-9]+")) {
				out.add(ctx);
				continue;
			}

			final List<Allele> alts = ctx.getAlternateAlleles();

			if (!AcidNucleics.isATGC(ctx.getReference())
					|| alts.stream().noneMatch(A -> AcidNucleics.isATGC(A))) {
				out.add(ctx);
				continue;
			}

			for (PopulationData pop : popHash.values()) {
				pop.reset(alts.size());
			}

			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			boolean failed_flag = false;
			for (int i = 0; i < alts.size(); ++i) {
				final Allele alt = alts.get(i);
				if (!AcidNucleics.isATGC(alt)) {
					// ignore
				} else {

					Response resp = null;
					if (debug)
						LOG.info(ctx.getContig() + ":" + ctx.getStart() + ":"
								+ ctx.getReference().getBaseString() + ":" + alt.getBaseString());
					for (int n_try = 0; n_try < this.max_retries; ++n_try) {
						final long wait_millisec = (this.max_seconds * 1000L)
								- (System.currentTimeMillis() - last_call);
						if (wait_millisec > 0L) {
							try {
								if (debug)
									LOG.info("wait " + wait_millisec + " ms");
								Thread.sleep(wait_millisec);
							} catch (InterruptedException ie) {
							}
						}
						resp = query(ctx, alt);
						last_call = System.currentTimeMillis();

						if (!resp.fail_flag)
							break;
						if (debug)
							LOG.info("trying " + (n_try + 2) + "/" + this.max_retries);
						resp = null;
					}
					if (resp == null) {
						failed_flag = true;
						break;
					} else if (resp.fail_flag || resp.root == null || !resp.root.isJsonObject()
							|| resp.root.getAsJsonObject().has("errors")) {
						if (resp.fail_flag) {
							failed_flag = true;
							break;
						}
					} else {
						JsonObject obj = resp.root.getAsJsonObject();
						if (obj.has("data") && obj.get("data").isJsonObject()) {
							obj = obj.get("data").getAsJsonObject();
							if (obj.has("variant") && obj.get("variant").isJsonObject()) {
								obj = obj.get("variant").getAsJsonObject();
								populationDetails(obj, i, popHash.get(GENERAL_POP));

								if (obj.has("populationDetails")) {
									populationDetails(obj.get("populationDetails"), i, popHash);
								}
							}
						}
					}
				}

			}
			if (failed_flag) {
				vcb.attribute(info_ukbbstatus.getID(), true);
			} else {
				for (PopulationData pop : popHash.values()) {
					build(pop, vcb);
				}
			}
			out.add(vcb.make());
		}
		return 0;
	}

	@Override
	protected Logger getLogger() {
		return LOG;
	}

	@Override
	protected int beforeVcf() {
		if (!this.debug && Runtime.version().feature() < 23) {
			LOG.error("I don't know why but this tool doesn't work with 'old' jdk, ok with 23");
			return -1;
		}

		return super.beforeVcf();
	}

	@Override
	protected void afterVcf() {
		/* autocloseable in jdk23 but not in jdk17 */
		CloserUtil.close(this.httpClient);
		super.afterVcf();
	}

	public static void main(final String[] args) {
		new VcfUkbiobank().instanceMainWithExit(args);
	}

}
