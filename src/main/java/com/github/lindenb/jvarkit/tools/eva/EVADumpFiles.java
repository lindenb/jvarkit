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
package com.github.lindenb.jvarkit.tools.eva;

import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;

import org.apache.http.HttpException;
import org.apache.http.HttpStatus;
import org.apache.http.StatusLine;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.net.HttpStatusException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import htsjdk.samtools.util.StringUtil;

/**
BEGIN_DOC

### Example
```
$ java -jar dist/jvarkit.jar evadumpfiles --taxons 'cfamiliaris 31' 2> /dev/null |  verticalize  

>>> 2
$1	studyId	PRJEB56315
$2	studyName	Idiopathic Epilepsy in the Dutch partridge dog
$3	assemblyAccession	GCA_000002285.2
$4	assemblyName	CanFam3.1
$5	taxonomyCommonName	Dog
$6	taxonomyScientificName	Canis lupus familiaris
$7	taxonomyId	9615
$8	taxonomyCode	cfamiliaris
$9	assemblyCode	31
$10	file	https://ftp.ebi.ac.uk/pub/databases/eva/PRJEB56315/variants_pass.vcf.gz
<<< 2

>>> 3
$1	studyId	PRJEB56211
$2	studyName	Performance of Variant Pathogenicity Prediction Methods
$3	assemblyAccession	GCA_000002285.2
$4	assemblyName	CanFam3.1
$5	taxonomyCommonName	Dog
$6	taxonomyScientificName	Canis lupus familiaris
$7	taxonomyId	9615
$8	taxonomyCode	cfamiliaris
$9	assemblyCode	31
$10	file	https://ftp.ebi.ac.uk/pub/databases/eva/PRJEB56211/VPP_Project_Variants.vcf.gz
<<< 3

>>> 4
$1	studyId	PRJEB51024
$2	studyName	Longevity of Cane Corso Italiano purebred dogs
$3	assemblyAccession	GCA_000002285.2
$4	assemblyName	CanFam3.1
$5	taxonomyCommonName	Dog
$6	taxonomyScientificName	Canis lupus familiaris
$7	taxonomyId	9615
$8	taxonomyCode	cfamiliaris
$9	assemblyCode	31
$10	file	https://ftp.ebi.ac.uk/pub/databases/eva/PRJEB51024/cane_corso_f.vcf.gz
```

END_DOC

*/

@Program(name="evadumpfiles",
keywords={"eva","ebi","snp","variant"},
description="Dump files locations from European Variation Archive",
modificationDate="20230314",
creationDate="20230314",
jvarkit_amalgamion = true
)
public class EVADumpFiles extends Launcher {
	private static final Logger LOG = Logger.of(EVADumpFiles.class);

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--taxons","--taxon"},description="define one or more taxon filters:syntax 'taxonCode assemblyCode' or 'taxonCode *' or 'taxonCode' .eg: 'cfamiliaris' , 'cfamiliaris * ', 'cfamiliaris 31'. Multiple are comma/semicolon separated. See output of https://www.ebi.ac.uk/eva/webservices/rest/v1/meta/species/list/ ." )
	private List<String> taxon_filter_str = new ArrayList<>();

	
	private CloseableHttpClient httpClient = null;
	private final JsonParser jsonparser = new JsonParser();
	private final List<AssemblyFilter> assembyFilters = new ArrayList<>();
	
	@SuppressWarnings("unused")
	private static class Assembly {
		String assemblyAccession;
	    String assemblyChain;
	    String assemblyVersion;
	    String assemblyName;
	    String assemblyCode;
	    
	    int taxonomyId;
	    String taxonomyCommonName;
	    String taxonomyScientificName;
	    String taxonomyCode;
	    String taxonomyEvaName;
	    
	    public String getTaxonomyAndAssembly() {
	    	return taxonomyCode+"_"+assemblyCode;
	    }
	    
	    @Override
	    public String toString() {
	    	return taxonomyCommonName;
	    	}
		}
	
	private static class Study {
		Assembly assembly;
		String studyId;
		String studyName;
		}
	
	private JsonObject fetchJSON(final String urlstr) throws IOException,HttpException {
		HttpGet httpGet = null;
		try {
			LOG.info(urlstr);
			httpGet = new HttpGet(urlstr);
			final CloseableHttpResponse httpResponse = httpClient.execute(httpGet);
			final StatusLine statusLine = httpResponse.getStatusLine();
			final int responseCode = statusLine.getStatusCode();
		  
			if(responseCode == HttpStatus.SC_NO_CONTENT || responseCode==HttpStatus.SC_NOT_FOUND /* yeah.. sometimes it happens.. */)
		 		{
				throw new HttpStatusException(statusLine);
		 		}
			else if(responseCode !=  HttpStatus.SC_OK)
			 	{
				throw new IOException("Remote error (http status: "+responseCode+") "+urlstr);
			 	}
			 else
				 {
			 	//response = new TeeInputStream( httpConnection.getInputStream(),System.err,false);
			 	try(InputStream response =httpResponse.getEntity().getContent()) {
					final JsonElement root = this.jsonparser.parse(new InputStreamReader(response));
					if(!root.isJsonObject()) {
						throw new IOException("not a json object: "+urlstr);
						}
					return root.getAsJsonObject();
					}
				 }
			}
		catch(final HttpException err) {
			LOG.error(err);
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
	
	private static class AssemblyFilter implements Predicate<Assembly> {
		final String taxonomyCode;
		final String assemblyCode;
		AssemblyFilter(String s) {
			final String[] tokens= s.split("\\s+");
			if(tokens.length>2) throw new IllegalArgumentException("expected one or two tokens in "+s);
			taxonomyCode  = tokens[0];
			if(tokens.length>1 && !tokens[1].equals("*") ) {
				assemblyCode = tokens[1];
				}
			else
				{
				assemblyCode = null;
				}
			}
		@Override
		public boolean test(Assembly t) {
			if(!t.taxonomyCode.equalsIgnoreCase(this.taxonomyCode)) return false;
			if(this.assemblyCode==null) return true;
			return t.assemblyCode.equalsIgnoreCase(this.assemblyCode);
			}
		}

	
	private List<Assembly> fetchAssemblies() throws IOException,HttpException {
		final String urlstr = "https://www.ebi.ac.uk/eva/webservices/rest/v1/meta/species/list/";/* WARNING it's v1 not v2 */
		final List<Assembly> L = new ArrayList<>();
		for(JsonElement e0 :fetchJSON(urlstr).get("response").getAsJsonArray()) {
			for(JsonElement e: e0.getAsJsonObject().get("result").getAsJsonArray()) {
				final JsonObject eo = e.getAsJsonObject();
				final Assembly as = new Assembly();
				as.assemblyAccession = eo.get("assemblyAccession").getAsString();
				as.assemblyChain = eo.get("assemblyChain").getAsString();
				as.assemblyVersion = eo.get("assemblyVersion").getAsString();
				as.assemblyName = eo.get("assemblyName").getAsString();
				as.assemblyCode = eo.get("assemblyCode").getAsString();
				as.taxonomyId = eo.get("taxonomyId").getAsInt();
				as.taxonomyCommonName = eo.get("taxonomyCommonName").getAsString();
				as.taxonomyScientificName = eo.get("taxonomyScientificName").getAsString();
				as.taxonomyCode = eo.get("taxonomyCode").getAsString();
				as.taxonomyEvaName = eo.get("taxonomyEvaName").getAsString();
				if(!this.assembyFilters.isEmpty() && this.assembyFilters.stream().noneMatch(AF->AF.test(as))) continue;
				L.add(as);
				}
			}
		return L;
		}
	private List<Study> fetchStudies() throws IOException,HttpException {
		List<Study> L = new ArrayList<>();
		for(Assembly as:fetchAssemblies()) {
			LOG.info("assembly :"+as);
			String urlstr = "https://www.ebi.ac.uk/eva/webservices/rest/v2/studies?species="+as.taxonomyCode+"&assembly="+as.assemblyCode;
			for(;;) {
				JsonObject root = null;
				
				try {
					root =  fetchJSON(urlstr).getAsJsonObject();
					}
				catch(HttpStatusException err) {
					continue;
					}
				for(JsonElement e:root.get("_embedded").getAsJsonObject().get("variantStudySummaryList").getAsJsonArray()) {
					final JsonObject eo = e.getAsJsonObject();
					if(eo.get("filesCount").getAsInt()==0) continue;
					Study study = new Study();
					study.assembly = as;
					study.studyId = eo.get("studyId").getAsString();
					study.studyName = eo.get("studyName").getAsString();
					L.add(study);
					}
				if(!root.has("_links")) break;
				JsonObject x = root.getAsJsonObject("_links");
				if(!x.has("next")) break;
				x = x.getAsJsonObject("next");
				if(!x.has("href")) break;
				urlstr = x.get("href").getAsString();
				}
			}
		return L;
		}
	
	private void fetchFiles(PrintWriter out) throws IOException,HttpException {
		final List<Study> studies =fetchStudies();
		LOG.info("N studies="+studies.size());
		for(Study study :studies) {
		final String urlstr="https://www.ebi.ac.uk/eva/webservices/rest/v1/studies/" + study.studyId +
				"/files?species="+ study.assembly.getTaxonomyAndAssembly();
		
		JsonArray array;
		try {
			array = fetchJSON(urlstr).get("response").getAsJsonArray();
			}
		catch(HttpException err) {
			LOG.error(err);
			continue;
			}
		for(JsonElement  e0:  array) {
			for(JsonElement  e1:  e0.getAsJsonObject().get("result").getAsJsonArray()) {
				final JsonObject eo = e1.getAsJsonObject();
				final String fileName = eo.get("fileName").getAsString();
				out.print(study.studyId);
				out.print("\t");
				out.print(study.studyName);
				out.print("\t");
				out.print(study.assembly.assemblyAccession);
				out.print("\t");
				out.print(study.assembly.assemblyName);
				out.print("\t");
				out.print(study.assembly.taxonomyCommonName);
				out.print("\t");
				out.print(study.assembly.taxonomyScientificName);
				out.print("\t");
				out.print(study.assembly.taxonomyId);
				out.print("\t");
				out.print(study.assembly.taxonomyCode);
				out.print("\t");
				out.print(study.assembly.assemblyCode);
				out.print("\t");
				out.println("https://ftp.ebi.ac.uk/pub/databases/eva/"+study.studyId+"/"+fileName);
				}
			}
		}
	}
	@Override
	public int doWork(final List<String> args) {
		if(!args.isEmpty()) {
			LOG.error("too many arguments");
			return -1;
			}
		try {
			this.taxon_filter_str.stream().
				flatMap(S->Arrays.stream(S.split("[,;]"))).
				filter(S->!StringUtil.isBlank(S)).
				map(S->new AssemblyFilter(S)).
				forEach(F->assembyFilters.add(F));
			
		
			this.httpClient = HttpClients.createSystem();//createDefault();
			try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				out.print("studyId");
				out.print("\t");
				out.print("studyName");
				out.print("\t");
				out.print("assemblyAccession");
				out.print("\t");
				out.print("assemblyName");
				out.print("\t");
				out.print("taxonomyCommonName");
				out.print("\t");
				out.print("taxonomyScientificName");
				out.print("\t");
				out.print("taxonomyId");
				out.print("\t");
				out.print("taxonomyCode");
				out.print("\t");
				out.print("assemblyCode");
				out.print("\t");
				out.println("file");
				fetchFiles(out);
				out.flush();
				}
			return 0;
		} catch(Throwable err) {
			LOG.error(err);
			return -1;
		}
		finally {
			if(this.httpClient!=null) try { httpClient.close();} catch(Throwable err2) {}
		}
	}
	
	
	public static void main(String[] args) {
		new EVADumpFiles().instanceMainWithExit(args);

	}
}
