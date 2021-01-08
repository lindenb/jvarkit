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
package com.github.lindenb.jvarkit.tools.spliceai;

import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

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

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonParser;
/**
BEGIN_DOC

# Example

```
java -jar dist/vcfspliceai.jar  src/test/resources/test_vcf01.vcf 

(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
(...)
1	866893	.	T	C	431	PASS	AA=t;AC=7;AF=0.7;AN=10;SPLICEAI=SAMD11|0.00|0.00|0.00|0.00|-13|24|-13|-48
1	870317	.	G	A	12	PASS	AC=11;AF=0.917;AN=12;SPLICEAI=SAMD11|0.00|0.00|0.00|0.00|2|17|16|-12
1	875770	.	A	G	338	PASS	AA=a;AC=8;AF=0.8;AN=10;SPLICEAI=SAMD11|0.00|0.00|0.01|0.00|-1|-45|-50|-46
1	903245	.	A	G	199	PASS	AA=a;AC=6;AF=0.6;AN=10;SPLICEAI=PLEKHN1|0.00|0.00|0.00|0.00|48|-37|-22|1
1	905130	.	ATG	A	487	PASS	AC=3;AF=0.5;AN=6;CIGAR=1M2D;IDREP=1;REFREP=2;RU=TG;SPLICEAI=PLEKHN1|0.00|0.00|0.00|0.00|-43|21|-33|-37
1	909238	.	G	C	229	PASS	AA=C;AC=8;AF=0.667;AN=12;SPLICEAI=PLEKHN1|0.00|0.01|0.00|0.00|-43|-50|39|-7
1	912049	.	T	C	400	PASS	AA=T;AC=5;AF=0.625;AN=8;SPLICEAI=PERM1|0.00|0.01|0.01|0.00|-28|-14|-27|-23
1	913889	.	G	A	372	PASS	AA=G;AC=5;AF=0.625;AN=8;SPLICEAI=PERM1|0.00|0.01|0.00|0.00|-46|9|2|-45
1	914333	.	C	G	556	PASS	AA=G;AC=5;AF=0.625;AN=8;SPLICEAI=PERM1|0.00|0.00|0.00|0.00|-3|27|-3|-38
1	914852	.	G	C	525	PASS	AA=C;AC=5;AF=0.625;AN=8;SPLICEAI=PERM1|0.00|0.00|0.00|0.00|22|-22|48|49
1	914940	.	T	C	488	PASS	AA=C;AC=5;AF=0.625;AN=8;SPLICEAI=PERM1|0.00|0.00|0.00|0.00|28|-30|-39|3
(...)
```

END_DOC
 */
@Program(name="vcfspliceai",
description="Annotate VCF with spiceai web service",
keywords={"vcf","splice","splicing","spliceai"},
creationDate="20201107",
modificationDate="20201107"
)
public class VcfSpliceAI  extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.build(VcfSpliceAI.class).make();
	
	private CloseableHttpClient httpClient = null;
	private final JsonParser jsonparser = new JsonParser();
	
	@Parameter(names={"--tag"},description="INFO tag")
	private String tag="SPLICEAI";
	@Parameter(names={"--base"},description="Base API")
	private String base="https://spliceailookup-api.broadinstitute.org/spliceai/";
	@Parameter(names={"--hg"},description="genome version. Must be 37 or 38. Otherwise, the dictionary is used to detect the version.")
	private int build=-1;
	@Parameter(names={"--distance"},description="For each variant, SpliceAI looks within a window (+/- 50bp by default) to see how the variant affects the probabilities of different positions being splice acceptors or donors. The distance specified here controls the size of this window. The maximum allowed value is 10,000bp")
	private int distance=50;
	@Parameter(names={"--raw"},description="Splicing changes corresponding to strengthening annotated splice sites and weakening unannotated splice sites are typically much less pathogenic than weakening annotated splice sites and strengthening unannotated splice sites. Selecting 'masked' (default) will hide the score for such splicing changes and show 0 instead. Selecting 'raw' will show all scores. SpliceAI developers recommend using 'raw' scores for alternative splicing analysis and 'masked' scores for variant interpretation.")
	private boolean raw = false;

	@Override
	protected Logger getLogger()
		{
		return LOG;
		}

	@Override
	protected int beforeVcf()
		{
		/** create http client */
		this.httpClient = HttpClients.createSystem();//createDefault();
		return super.beforeVcf();
		}
	
	@Override
	protected void afterVcf() {
		CloserUtil.close(this.httpClient);
		this.httpClient=null;
		super.afterVcf();
		}
	
	private static class SpliceResponse {
		String errorMsg;
		String scores ; 
	}
	
	private SpliceResponse callApi(final String contig,int start,Allele ref,Allele alt) {
		
		InputStream response =null;
		HttpGet httpGet = null;
		final SpliceResponse apiReturn = new SpliceResponse();
		
		if(alt.equals(Allele.SPAN_DEL) || alt.isSymbolic() || ref.isSymbolic()) {
			apiReturn.scores=".";
			return apiReturn;
			}
		try {
			final String urlstr = new StringBuilder(this.base).append("?hg=").
					append(this.build).
					append("&distance=").append(this.distance).
					append("&mask=").append(this.raw?0:1).
					append("&variant=").
					append(contig).append("-").
					append(start).append("-").
					append(ref.getDisplayString()).append("-").
					append(alt.getDisplayString()).toString();
			
			httpGet = new HttpGet(urlstr);
			final CloseableHttpResponse httpResponse = httpClient.execute(httpGet);
			final int responseCode = httpResponse.getStatusLine().getStatusCode();
		  
			 if(responseCode != 200)
			 	{
				apiReturn.errorMsg="HttpError"+responseCode;
			 	}
			 else
				 {
			 	//response = new TeeInputStream( httpConnection.getInputStream(),System.err,false);
			 	response =httpResponse.getEntity().getContent();
				final JsonElement root = this.jsonparser.parse(new InputStreamReader(response));
				if(root.getAsJsonObject().has("error")) {
					apiReturn.errorMsg = root.getAsJsonObject().get("error").getAsString();
					return apiReturn;
					}
				if(root.getAsJsonObject().has("scores")) {
					JsonArray array= root.getAsJsonObject().getAsJsonArray("scores");
					if(array.size()!=1) {
						LOG.error("illegal state "+urlstr);
						apiReturn.errorMsg="more than one score for "+urlstr;
						return apiReturn;
						}
					apiReturn.scores = array.get(0).getAsString();
					return apiReturn;
					}
				}
		
			}
		catch(final Throwable err) {
			LOG.error(err);
			apiReturn.errorMsg=String.valueOf(err.getMessage());
			}
		 finally
			{
			CloserUtil.close(response);
			if(httpGet!=null) httpGet.releaseConnection();
			}
		return apiReturn;
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin,
			VariantContextWriter out)
		{
		final VCFHeader header = iterin.getHeader();
		if( build==-1) {
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
		
			if(SequenceDictionaryUtils.isGRCh37(dict)) {
				build= 37;
				}
			else if(SequenceDictionaryUtils.isGRCh38(dict)) {
				build= 38;
				}
			}
		if(!(build==37 || build==38))
			{
			LOG.error("build hg"+build+" is not either grch37 or grch38.");
			return -1;
			}
		final VCFInfoHeaderLine info = new VCFInfoHeaderLine(this.tag, VCFHeaderLineCount.A, VCFHeaderLineType.String,"annotation from "+this.base);
		header.addMetaDataLine(info);
		JVarkitVersion.getInstance().addMetaData(this, header);
		final List<SpliceResponse> responses= new ArrayList<>();
		out.writeHeader(header);
		while(iterin.hasNext()) {
			final VariantContext ctx= iterin.next();
			if(!ctx.isVariant()) {
				out.add(ctx);
				continue;
				}
			responses.clear();
			for(final Allele alt:ctx.getAlternateAlleles()) {
				responses.add(callApi(ctx.getContig(), ctx.getStart(), ctx.getReference(), alt));
				}
			if(responses.stream().allMatch(R->R.scores!=null && R.scores.equals(".")) || 
				responses.stream().allMatch(R->R.errorMsg!=null)) {
				out.add(ctx);
				continue;
				}
			out.add(new VariantContextBuilder(ctx).attribute(info.getID(),
					responses.stream().map(R->R.errorMsg!=null?".":R.scores).collect(Collectors.toList())
					).make());
			}
		return 0;
		}
	
	public static void main(final String[] args) {
		new VcfSpliceAI().instanceMainWithExit(args);
		}
	
	}
