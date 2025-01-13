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
package com.github.lindenb.jvarkit.tools.gtex;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;

import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;

/**
BEGIN_DOC
## Motivation

export gtex eqtl data from a list of RS using the GTEX API

## Example
```
$ cat mylist.of.rs.txt | java -jar dist/jvarkit.jar gtexrs2qtl  | head | column -t

method            chromosome  datasetId  gencodeId           geneSymbol  geneSymbolUpper  nes        pValue       pos       snpId      tissueSiteDetailId                   variantId               phenotypeId
singleTissueEqtl  chr20       gtex_v8    ENSG00000088298.12  EDEM2       EDEM2            -0.292217  1.86762e-15  35045523  rs6088690  Heart_Left_Ventricle                 chr20_35045523_A_G_b38  .
singleTissueEqtl  chr20       gtex_v8    ENSG00000101000.5   PROCR       PROCR            -0.341834  9.7875e-10   35045523  rs6088690  Heart_Left_Ventricle                 chr20_35045523_A_G_b38  .
singleTissueEqtl  chr20       gtex_v8    ENSG00000088298.12  EDEM2       EDEM2            -0.360261  6.66882e-07  35045523  rs6088690  Brain_Putamen_basal_ganglia          chr20_35045523_A_G_b38  .
singleTissueEqtl  chr20       gtex_v8    ENSG00000101460.12  MAP1LC3A    MAP1LC3A         0.207036   3.12021e-10  35045523  rs6088690  Thyroid                              chr20_35045523_A_G_b38  .
singleTissueEqtl  chr20       gtex_v8    ENSG00000078814.15  MYH7B       MYH7B            0.217765   3.34276e-17  35045523  rs6088690  Thyroid                              chr20_35045523_A_G_b38  .
singleTissueEqtl  chr20       gtex_v8    ENSG00000100991.11  TRPC4AP     TRPC4AP          -0.174607  4.81976e-12  35045523  rs6088690  Thyroid                              chr20_35045523_A_G_b38  .
singleTissueEqtl  chr20       gtex_v8    ENSG00000088298.12  EDEM2       EDEM2            -0.265171  2.8512e-15   35045523  rs6088690  Thyroid                              chr20_35045523_A_G_b38  .
singleTissueEqtl  chr20       gtex_v8    ENSG00000101000.5   PROCR       PROCR            -0.181498  9.80399e-07  35045523  rs6088690  Thyroid                              chr20_35045523_A_G_b38  .
singleTissueEqtl  chr20       gtex_v8    ENSG00000078814.15  MYH7B       MYH7B            0.211821   3.58767e-08  35045523  rs6088690  Esophagus_Gastroesophageal_Junction  chr20_35045523_A_G_b38  .
```

END_DOC
*/
@Program(name="gtexrs2qtl",
description="extract gtex eqtl data from a list of RS",
keywords={"gtex","rs","eqtl","sqtl"},
creationDate="20230215",
modificationDate="20240225",
jvarkit_amalgamion = true
)
public class GtexRsToQTL extends Launcher {
private static final Logger LOG=Logger.build(GtexRsToQTL.class).make(); 
@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
private Path outputFile = null;
@Parameter(names={"--base-api"},description="GTEX base API")
private String gtex_base_api = "https://gtexportal.org/api/v2";

@Parameter(names={"--debug-json"},hidden = true)
private boolean debug_json = false;

private static final String[] TISSUE_COLS = new String[] {
	"chromosome",
	"datasetId",
	"gencodeId",
	"geneSymbol",
	"geneSymbolUpper",
	"nes",
	"pValue",
	"pos",
	"snpId",
	"tissueSiteDetailId",
	"variantId",
	"phenotypeId"
	};

private static final class Params {
	CloseableHttpClient httpClient;
	final JsonParser jsonparser = new JsonParser();
	PrintWriter out;
	}

private void test_paginf(JsonElement e) {
	if(!e.isJsonObject()) return;
	JsonObject o = e.getAsJsonObject();
	if(!o.has("paging_info")) return;
	e= o.get("paging_info");
	if(!e.isJsonObject()) return;
	o = e.getAsJsonObject();
	if(!o.has("numberOfPages")) return;
	if(o.get("numberOfPages").getAsInt()>1) throw new IllegalStateException("fix this numpage>1");
	}

private JsonElement call(
		final Params params,
		final String urlstr) throws IOException {
	LOG.info(urlstr);
	HttpGet httpGet=null;
	try {
		httpGet = new HttpGet(urlstr);
		final CloseableHttpResponse httpResponse = params.httpClient.execute(httpGet);
		final int responseCode = httpResponse.getStatusLine().getStatusCode();
	  
		 if(responseCode != 200)
		 	{
			throw new IOException("HttpError:"+responseCode+" for "+urlstr);
		 	}
		 else
			 {
		 	//response = new TeeInputStream( httpConnection.getInputStream(),System.err,false);
			 try(InputStream response =httpResponse.getEntity().getContent()) {
				final String json =  IOUtil.slurp(response);
				if(debug_json ) LOG.warn(json);
				return params.jsonparser.parse(json);
			 	}
			}
	
		}
	catch(final IOException err) {
		throw err;
		}
	catch(final Throwable err) {
		throw new IOException(err);
		}
	 finally
		{
		if(httpGet!=null) httpGet.releaseConnection();
		}
	}

private void tissue(
		final Params params,
		final String method,
		final String variantId
		)  throws IOException {
	final String url= this.gtex_base_api + "/association/"+ method +
			"?format=json&variantId=" + variantId + "&datasetId=gtex_v8";
	final PrintWriter out = params.out;
	final JsonElement e  = call(params,url);
	test_paginf(e);
	final JsonObject object1 = e.getAsJsonObject();
	if(!object1.has("data")) {
		LOG.warning("no data in "+e);
		return;
		}
	final JsonArray array = object1.get("data").getAsJsonArray();
	for(JsonElement e2:array) {
		final JsonObject object2 = e2.getAsJsonObject();
		out.print(method);
		for(int i=0;i< TISSUE_COLS.length;i++) {
			final String col= TISSUE_COLS[i];
			out.print("\t");
			if(!object2.has(col))
				{
				out.print(".");
				}
			else
				{
				out.print(object2.get(col).getAsString());
				}
			}
		out.println();
		object2.entrySet().stream().map(KV->KV.getKey()).
			filter(S->Arrays.stream(TISSUE_COLS).noneMatch(C->C.equals(S))).
			forEach(S->LOG.warn("Col "+S+" missing for "+method));
		}
	}

private void oneVariantId(
		final Params params, 
		final String variantId
		) throws IOException {
	tissue(params,"singleTissueEqtl",variantId);
	tissue(params,"singleTissueSqtl",variantId);
	}

private void oneRS(
		final Params params, 
		final String rs
		) throws IOException {
	JsonElement e1 = call(params, this.gtex_base_api + "/dataset/variant?format=json&snpId="+ rs + "&datasetId=gtex_v8");
	test_paginf(e1);
	JsonObject o1 = e1.getAsJsonObject();
	if(!o1.has("data")) return;
	
	JsonArray array = o1.get("data").getAsJsonArray();
	if(array.size()==0) return;
	for(JsonElement e2: array) {
		JsonObject o2 = e2.getAsJsonObject();
		if(!o2.has("variantId")) continue;
		final String variantId = o2.get("variantId").getAsString();
		if(StringUtils.isBlank(variantId)) continue;
		oneVariantId(params,variantId);
		}
	}

@Override
public int doWork(final List<String> args) {
	final Params params = new Params();
	try {
		try(PrintWriter pw= super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
			params.out = pw;
			params.out.println("method\t"+String.join("\t",TISSUE_COLS));
			try(CloseableHttpClient httpClient = HttpClients.createSystem()) {
				params.httpClient = httpClient;
				if(args.isEmpty()) {
					try(BufferedReader r= IOUtils.openStreamForBufferedReader(stdin())) {
						String line;
						while((line=r.readLine())!=null) {
							if(StringUtil.isBlank(line) || line.startsWith("#")) continue;
							oneRS(params,line.trim());
							}
						}
					}
				else
					{
					for(String rs: args) {
						oneRS(params,rs);
						}
					}
				}
			params.out.flush();
			}
		return 0;
		}
	catch(final Throwable err) {
		LOG.error(err);
		return -1;
		}
	}	
	
public static void main(final String[] args) {
	new GtexRsToQTL().instanceMainWithExit(args);
	}
}
