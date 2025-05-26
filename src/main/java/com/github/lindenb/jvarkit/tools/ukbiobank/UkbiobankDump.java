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
package com.github.lindenb.jvarkit.tools.ukbiobank;

import java.io.IOException;
import java.io.PrintWriter;
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
import java.nio.file.Path;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Locatable;
/**
BEGIN_DOC

# WARNING

I cannot explain why this software doesn't run with jdk17 (I always get a http error 403) but runs with jdk23...



```
$ java -jar jvarkit.jar ukbbdump -R ref.dict --region "chr3:38000000-39000000" 
Chrom	Pos	rsID	Ref	Alt	nHomozygotes	HGVSp	maxImpact	maxConsequence	alleleCount	alleleNum	alleleFreq	geneSymbol
chr3	38000013	.	A	G	0	NM_001370264.1:c.672+574A>G,NM_001370265.1:c.573+574A>G,NM_001385038.1:c.1182+574A>G,NM_001385039.1:c.1182+574A>G,NM_015873.4:c.1182+574A>G	LOWEST	intron_variant	13	981062	0.0000132509464233657	VILL
chr3	38000014	.	G	A	0	NM_001370264.1:c.672+575G>A,NM_001370265.1:c.573+575G>A,NM_001385038.1:c.1182+575G>A,NM_001385039.1:c.1182+575G>A,NM_015873.4:c.1182+575G>A	LOWEST	intron_variant	3	981092	0.0000030578172077644096	VILL
chr3	38000015	.	G	A	0	NM_001370264.1:c.672+576G>A,NM_001370265.1:c.573+576G>A,NM_001385038.1:c.1182+576G>A,NM_001385039.1:c.1182+576G>A,NM_015873.4:c.1182+576G>A	LOWEST	intron_variant	3	981092	0.0000030578172077644096	VILL
chr3	38000018	rs1319217247	G	A	0	NM_001370264.1:c.672+579G>A,NM_001370265.1:c.573+579G>A,NM_001385038.1:c.1182+579G>A,NM_001385039.1:c.1182+579G>A,NM_015873.4:c.1182+579G>A	LOWEST	intron_variant	2	981088	0.0000020385531165400047	VILL
chr3	38000021	rs6806791	G	A	32	NM_001370264.1:c.672+582G>A,NM_001370265.1:c.573+582G>A,NM_001385038.1:c.1182+582G>A,NM_001385039.1:c.1182+582G>A,NM_015873.4:c.1182+582G>A	LOWEST	intron_variant	1342	981080	0.0013678802951848984	VILL
chr3	38000022	rs937502348	T	C	0	NM_001370264.1:c.672+583T>C,NM_001370265.1:c.573+583T>C,NM_001385038.1:c.1182+583T>C,NM_001385039.1:c.1182+583T>C,NM_015873.4:c.1182+583T>C	LOWEST	intron_variant	4	981070	0.0000040771810370309966	VILL
chr3	38000027	.	T	C	0	NM_001370264.1:c.672+588T>C,NM_001370265.1:c.573+588T>C,NM_001385038.1:c.1182+588T>C,NM_001385039.1:c.1182+588T>C,NM_015873.4:c.1182+588T>C	LOWEST	intron_variant	1	981072	0.0000010192931813363342	VILL
chr3	38000033	.	C	G	0	NM_001370264.1:c.672+594C>G,NM_001370265.1:c.573+594C>G,NM_001385038.1:c.1182+594C>G,NM_001385039.1:c.1182+594C>G,NM_015873.4:c.1182+594C>G	LOWEST	intron_variant	3	981056	0.000003057929414834627	VILL
```
##

END_DOC

 */
@Program(
		name="ukbbdump",
		description="Dump data https://afb.ukbiobank.ac.uk/ ukbiobank server. The server might not like too many requests. Use a your own risk. Doesn't work with jdk17 (?!)",
		keywords={"ukbiobank"},
		creationDate="20250425",
		modificationDate="20250425",
		generate_doc = true,
		jvarkit_amalgamion = true
		)
public class UkbiobankDump extends Launcher {
	private static final Logger LOG = Logger.of(UkbiobankDump.class);
	
	private final JsonParser jsonparser = new JsonParser();
	@Parameter(names = { "-o" }, description = OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile =null;

	@Parameter(names = { "-R","--dict" }, description = DICTIONARY_SOURCE, required=true)
	private Path dictPath=null;
	@Parameter(names = { "--region" }, description = "limit to that region. "+IntervalParser.OPT_DESC )
	private String regionStr;
	@Parameter(names = { "-w" }, description = "window size",hidden = true)
	private int window_size = 100;

	
	@Parameter(names = { "--api" }, description = "API base url")
	private String api_url = "https://afb.ukbiobank.ac.uk/api";
	@Parameter(names = { "--seconds" }, description = "wait 'x' seconds between each call to the API")
	private int max_seconds = 10;
	@Parameter(names = { "--retry" }, description = "max-retry")
	private int max_retries = 10;
	@Parameter(names = { "--debug" }, description = "debug")
	private boolean debug = false;
	
	private long last_call=0L;
	private static final String[] HEADERS=new String[] {
		"Chrom","Pos","rsID","Ref","Alt","nHomozygotes","HGVSp","maxImpact","maxConsequence","alleleCount","alleleNum","alleleFreq","geneSymbol"
		};
	private static class Response {
		boolean fail_flag=false;
		JsonElement root;
		@Override
		public String toString() {
			return "fail:"+fail_flag+" root:"+root;
			}
	}
	
	private String getContig(final String ctg) {
		return ctg.startsWith("chr") ? ctg : "chr" + ctg;
	}



	private Response query(final HttpClient httpClient,String contig,int start,int end) {
		final Response apiReturn = new Response();
		
		final String query0 = "{\"query\":\"\nquery Query($regionStr: String!) {\n  region (regionStr: $regionStr) {\n    regionStr\n    variants {\n        rsID\n        Chrom\n        Pos\n        Ref\n        Alt\n        nHomozygotes\n        HGVSp\n        maxImpact\n        maxConsequence\n        alleleCount\n        alleleNum\n        alleleFreq\n        nHomozygotes\n        geneSymbol\n    }\n  }\n}\n\",\"variables\":{\"regionStr\":\""+contig+"-"+start+"-"+end+"\"}}";
		final String query = query0.replace("\n", "\\n");
		if (debug)
			LOG.info("sending "+query);
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
				hb.header("Referer", "https://afb.ukbiobank.ac.uk/region/" + contig + "-" + start+"-"+end);
				hb.uri(URI.create(this.api_url));
				hb.POST(HttpRequest.BodyPublishers.ofString(query, StandardCharsets.UTF_8));
				final HttpRequest request = hb.build();
				final HttpResponse<String> resp = httpClient.send(request,
						HttpResponse.BodyHandlers.ofString());

				if (resp.statusCode() == 200) {
					final String jsonStr = resp.body();
					if (jsonStr.startsWith("<")) {
						LOG.info("looks like html ");
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

	private int scan(HttpClient httpClient ,PrintWriter w, Locatable loc) {
		final String contig = getContig(loc.getContig());
		if (!contig.matches("(chr)?[XY0-9]+")) {
			return 0;
			}
		int x=loc.getStart();
		while(x<=loc.getEnd()) {
			final int x2 = Math.min(loc.getEnd(),x+this.window_size);
			
			Response resp = null;
			if (debug) LOG.info(contig+":"+x+"-"+x2);
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
				resp = query(httpClient, contig,x,x2);
				last_call = System.currentTimeMillis();

				if (!(resp.fail_flag || resp.root == null || !resp.root.isJsonObject()
						|| resp.root.getAsJsonObject().has("errors"))) {
					if (debug) LOG.info("Looks OK");
					break;
					}
				else
					{
					if (debug) LOG.info("Looks  BAD"+resp);
					}
				if (debug) LOG.info("trying " + (n_try + 2) + "/" + this.max_retries);
				resp = null;
				}
			if (resp == null || resp.fail_flag || resp.root == null || !resp.root.isJsonObject()
					|| resp.root.getAsJsonObject().has("errors")) {
				LOG.error("Error "+loc);
				return -1;
				} 
			JsonObject obj = resp.root.getAsJsonObject();
			if (obj.has("data") && obj.get("data").isJsonObject()) {
				obj = obj.get("data").getAsJsonObject();
				if (obj.has("region") && obj.get("region").isJsonObject()) {
					obj = obj.get("region").getAsJsonObject();
					if (obj.has("variants") && obj.get("variants").isJsonArray()) {
						JsonArray array = obj.get("variants").getAsJsonArray();
						for(int i=0;i< array.size();i++) {
							obj = array.get(i).getAsJsonObject();
							for(int j=0;j< HEADERS.length;j++) {
								if(j>0) w.print("\t");
								if(obj.has(HEADERS[j])) {
									w.print(obj.get(HEADERS[j]).getAsString());
									}
								}
							w.println();
						}
						
						
					}
					
					
				}
			}	
			
			x=x2+1;
		}
	return 0;
	}
	
	@Override
	public int doWork(List<String> args) {
		if(!args.isEmpty()) {
			LOG.error("no argument allowed");
			return -1;
			}
		
		
		if (!this.debug && Runtime.version().feature() < 23) {
			LOG.error("I don't know why but this tool doesn't work with 'old' jdk, ok with 23");
			return -1;
		}
		
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.dictPath);
		if (!SequenceDictionaryUtils.isGRCh38(dict)) {
			LOG.error("sequence dictionary doesn't look like hg38 ?!!:"
					+ SequenceDictionaryUtils.getBuildName(dict).orElse("undefined"));
			return -1;
			}
		HttpClient httpClient=null;
		try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
			httpClient   = HttpClient.newHttpClient();
			out.println(String.join("\t", HEADERS));
			
			
			if(!StringUtils.isBlank(this.regionStr)) {
				return scan(httpClient, out, new IntervalParser(dict).apply(this.regionStr).get());
				}
			else
				{
				for(SAMSequenceRecord ssr:dict.getSequences()) {
					if(scan(httpClient, out, ssr)!=0) return -1;
					}
				}
			out.flush();
			return 0;
			}
		catch(IOException err) {
			LOG.error(err);
			return -1;
			}
		finally {
			CloserUtil.close(httpClient);
			}
		}


	public static void main(final String[] args) {
		new UkbiobankDump().instanceMainWithExit(args);
	}

}
