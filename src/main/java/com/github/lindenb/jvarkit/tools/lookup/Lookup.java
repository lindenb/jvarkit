package com.github.lindenb.jvarkit.tools.lookup;

import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.http.HttpRequest;
import org.apache.http.client.HttpClient;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.methods.HttpPost;
import org.apache.http.client.methods.HttpRequestBase;
import org.apache.http.entity.ByteArrayEntity;
import org.apache.http.entity.ContentType;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.w3c.dom.Document;

import com.github.lindenb.jvarkit.io.TeeInputStream;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.tools.ensemblreg.VcfEnsemblReg;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;

public class Lookup extends Launcher {
	private static final int MAX_QUERY_REGION_SIZE=4_999_000;
	private static final Logger LOG = Logger.of(Lookup.class);
	private boolean teeResponse =false;
	private long lastMillisec = 0L;
	public Lookup() {
		}
	
	private abstract class JSONWrapper {
		JsonObject object;
		JSONWrapper(JsonObject object) {
			this.object= object;
			}
		JsonElement get(final String k) {
			return Objects.requireNonNull(this.object.get(k),"cannot find "+k+" in "+toString());
			}
		
		String getString(String key) {
			return get(key).getAsString();
			}
		int getInt(String key) {
			return get(key).getAsInt();
			}
		double getDouble(String key) {
			return get(key).getAsDouble();
			}
		@Override
		public String toString() {
			return object.toString();
			}
		}
	private class EnsemblGene extends JSONWrapper implements Locatable {
		EnsemblGene(JsonObject object) {
			super(object);
			}
		String getGeneId() {
			return getString("gene_id");
			}
		String getGeneName() {
			return getString("external_name");
			}
		String getType() {
			return getString("feature_type");
			}
		String getDescription() {
			return getString("description");
			}
		@Override
		public int hashCode() {
			return getGeneId().hashCode();
			}
		@Override
		public boolean equals(Object obj) {
			if(obj==this) return true;
			return (obj instanceof EnsemblGene) && getGeneId().equals(EnsemblGene.class.cast(obj).getGeneId());
			}
		int getVersion() {
			return getInt("version");
		}
		@Override
		public String getContig() {
			return getString("seq_region_name");
			}
		@Override
		public int getStart() {
			return getInt("start");
			}
		@Override
		public int getEnd() {
			return getInt("end");
			}
		}
	
	private class OpenTargetDisease extends JSONWrapper{
		OpenTargetDisease(JsonObject object) {
			super(object);
			}
		String getId() {
			return get("disease").getAsJsonObject().get("id").getAsString();
			}
		String getName() {
			return get("disease").getAsJsonObject().get("name").getAsString();
			}
		double getScore() {
			return get("score").getAsDouble();
			}
		@Override
		public int hashCode() {
			return getId().hashCode();
			}
		@Override
		public boolean equals(Object obj) {
			if(obj==this) return true;
			return (obj instanceof OpenTargetDisease) && getId().equals(OpenTargetDisease.class.cast(obj).getId());
			}
		
		}
	private class MolecularTrait extends JSONWrapper{
		MolecularTrait(JsonObject object) {
			super(object);
			}
		public double getPValue() {
			return getDouble("pvalue");
			}
		public String getStudyId() {
			return getString("study_id");
			}
		public String getDataSetId() {
			return getString("dataset_id");
			}
		}
	private Set<String> keySet(final JsonObject o) {
		return o.entrySet().stream().map(KV->KV.getKey()).collect(Collectors.toSet());
		}
	
	private InputStream tee(InputStream in) {
		return (this.teeResponse ?
			new TeeInputStream(in,stderr(),false):
				in
		 	);
		}
	/** send a pool of variants to VEP, returns the DOM document */
	private JsonElement callRequestToJson(final HttpClient httpClient,final HttpRequestBase httpRequest ) throws IOException
		{
		LOG.info(httpRequest.toString());
		httpRequest.setHeader("Accept",ContentType.APPLICATION_JSON.getMimeType());
		if ( this.lastMillisec!=-1L && this.lastMillisec+ 5000<  System.currentTimeMillis())
	    	{
	    	LOG.debug("waiting");
	    	try {Thread.sleep(1000);} catch(Exception err){}
	    	}
		try(CloseableHttpResponse httpResponse = (CloseableHttpResponse)httpClient.execute(httpRequest)) {
		final int responseCode = httpResponse.getStatusLine().getStatusCode();
		 if(responseCode != 200)
		 	{
			throw new RuntimeIOException("Response code was not 200. Detected response was "+responseCode+" for "+httpRequest);
		 	}
		 
		 //response = new TeeInputStream( httpConnection.getInputStream(),System.err,false);
			try(InputStream in =tee(httpResponse.getEntity().getContent())) {
				try(Reader r = new InputStreamReader(in, "UTF-8")) {
					JsonParser parser = new JsonParser();
					JsonElement e= parser.parse(r);
					return e;
					}
				}
			}
			finally {
				this.lastMillisec = System.currentTimeMillis(); 
			}
		
		}
	
	private JsonElement callOpenTargets(final HttpClient httpClient,JsonObject query) throws IOException {
		final String postBody = query.toString();
		final HttpPost httpPost = new HttpPost("https://api.platform.opentargets.org/api/v4/graphql");
		httpPost.setHeader("Content-Type",ContentType.APPLICATION_JSON.getMimeType());
		httpPost.setHeader("Accept",ContentType.APPLICATION_JSON.getMimeType());
		httpPost.setEntity(new ByteArrayEntity(postBody.getBytes(), ContentType.APPLICATION_JSON));
		return callRequestToJson(httpClient,httpPost);
		}
		 

	
	
	private Set<EnsemblGene> getEnsemblGenes(final HttpClient httpClient,final Locatable loc) throws IOException {
		final Set<EnsemblGene> genes= new HashSet<>();
		int x1=loc.getStart();
		while(x1 < loc.getEnd()) {
			int x2 = Math.min(x1+ MAX_QUERY_REGION_SIZE,loc.getEnd());
			final HttpGet httpGet = new HttpGet("https://rest.ensembl.org/overlap/region/human/"
					+ loc.getContig()+":"+x1+"-"+x2+"?feature=gene;content-type=application/json");
			final JsonArray array = this.callRequestToJson(httpClient, httpGet).getAsJsonArray();
			for(int i=0;i< array.size();i++) {
				genes.add(new EnsemblGene(array.get(i).getAsJsonObject()));
				}
			x1=x2;
			}
		return genes;
		}
	private Set<OpenTargetDisease> callOpenTargetsForGene(final HttpClient httpClient,EnsemblGene ensGene) throws IOException {
		final Set<OpenTargetDisease> diseases = new HashSet<>();
		JsonObject query = new JsonObject();
		query.addProperty("query", "query GeneDiseases($ensemblId: String!) "
				+ "{ target(ensemblId: $ensemblId) { id approvedSymbol associatedDiseases("
				+ "page: { index: 0, size: 50 }) "
				+ "{ rows { disease { id name dbXRefs } score datasourceScores { id score } } count } } }");
		final JsonObject o2 = new JsonObject();
		o2.addProperty("ensemblId", ensGene.getGeneId());
		query.add("variables",o2);
		JsonObject resp= callOpenTargets(httpClient, query).getAsJsonObject();
		//System.err.println(resp);
		JsonObject data = Objects.requireNonNull(resp.get("data")).getAsJsonObject();
		JsonObject target = Objects.requireNonNull(data.get("target")).getAsJsonObject();
		//System.err.println(target);
		JsonArray associatedDiseases = Objects.requireNonNull(target.get("associatedDiseases")).getAsJsonObject().get("rows").getAsJsonArray();
		for(JsonElement associatedDisease :associatedDiseases) {
				JsonObject disease= associatedDisease.getAsJsonObject();
				System.err.println(disease);
				diseases.add(new  OpenTargetDisease(disease));
				}
		return diseases;
		}
	
	private void getEQTLForGene(final HttpClient httpClient,final EnsemblGene gene) throws IOException {
		 HttpGet httpGet = new HttpGet("https://www.ebi.ac.uk/eqtl/api/v3/associations?gene_id="+gene.getGeneId());
		final JsonArray a = this.callRequestToJson(httpClient, httpGet).getAsJsonArray();
		Set<String> studies = new HashSet<String>();
		for(int i=0;i< a.size();++i) {
			MolecularTrait trait = new MolecularTrait(a.get(i).getAsJsonObject());
			
			studies.add(trait.getStudyId());
			}
		for(String study_id: studies) {
			httpGet = new HttpGet("https://www.ebi.ac.uk/eqtl/api/v3/studies/"+study_id);
			System.err.println("STUDY"+this.callRequestToJson(httpClient, httpGet));
			}
		}
	
	private void getGTexEQTLForGene(final HttpClient httpClient,final EnsemblGene gene) throws IOException {
		int page =0;
		for(;;) {
			HttpGet httpGet = new HttpGet("https://gtexportal.org/api/v2/association/singleTissueEqtl?page="+page+"&itemsPerPage=250&gencodeId="+gene.getGeneId()+"."+gene.getVersion());
			final JsonObject root = this.callRequestToJson(httpClient, httpGet).getAsJsonObject();
			JsonArray data = root.get("data").getAsJsonArray();
			if(data.size()==0) break;
			for(int i=0;i< data.size();i++) {
				JsonObject qtl = data.get(i).getAsJsonObject();
				System.err.println("GETX QTL "+qtl);
				}
			//int numberOfPage = root.get("paging_info").getAsJsonObject().get("totalNumberOfItems").getAsInt();
			++page;
			}
		
		}
	
	private void getPhewebForGene(final HttpClient httpClient,final EnsemblGene gene) throws IOException {
		try {
			final HttpGet httpGet = new HttpGet("https://clsa-pheweb.cerc-genomic-medicine.ca/api/v1/gene/"+gene.getGeneId());
			final JsonElement a = this.callRequestToJson(httpClient, httpGet);;
			System.err.println(a);
			}
		catch(Throwable err) {
			err.printStackTrace();
			}
		}
	
	@Override
	public int doWork(List<String> args) {
		try {
			try( CloseableHttpClient httpClient = HttpClients.createSystem(); ) {
				for(EnsemblGene ensGene: getEnsemblGenes(httpClient,new Interval("chr3",38541867,38567253)))  {
					if(ensGene.getGeneName().equals("SCN5A")==false) continue;
					System.err.println(ensGene);
					callOpenTargetsForGene(httpClient,ensGene);
					getEQTLForGene(httpClient,ensGene);
					getPhewebForGene(httpClient,ensGene);
					getGTexEQTLForGene(httpClient,ensGene);
					}
				}
			return 0;
			}
		catch(Throwable err) {
			err.printStackTrace();
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new Lookup().instanceMainWithExit(args);
	}
}
