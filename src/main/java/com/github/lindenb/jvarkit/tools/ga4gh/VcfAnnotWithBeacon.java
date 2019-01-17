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
package com.github.lindenb.jvarkit.tools.ga4gh;

import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.security.cert.CertificateException;
import java.security.cert.X509Certificate;
import java.util.HashSet;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import javax.net.ssl.SSLContext;

import org.apache.http.HttpHost;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.entity.ContentType;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.TeeInputStream;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import htsjdk.variant.vcf.VCFIterator;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.sleepycat.bind.tuple.StringBinding;
import com.sleepycat.bind.tuple.TupleBinding;
import com.sleepycat.bind.tuple.TupleInput;
import com.sleepycat.bind.tuple.TupleOutput;
import com.sleepycat.je.Database;
import com.sleepycat.je.DatabaseConfig;
import com.sleepycat.je.DatabaseEntry;
import com.sleepycat.je.Environment;
import com.sleepycat.je.EnvironmentConfig;
import com.sleepycat.je.LockMode;
import com.sleepycat.je.OperationStatus;
import com.sleepycat.je.Transaction;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
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
/* 
 BEGIN_DOC
 
 
 END_DOC
 */
@Program(
		name="vcfannotwithbeacon",
		description="Annotate a VCF with ga4gh beacon",
		keywords={"ga4gh","beacon","vcf","annotation"}
		)
public class VcfAnnotWithBeacon extends Launcher {
	private static final Logger LOG=Logger.build(VcfAnnotWithBeacon.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-B","--bdb"},description="Optional BerkeleyDB directory to store result. Avoid to make the same calls to beacon")
	private File bdbDir = null;
	@Parameter(names={"--build"},description="genome build")
	private String genomeBuild = "HG19";
	@Parameter(names={"--tag","-T"},description="INFO TAG")
	private String infoTag = "BEACON";
	@Parameter(names={"--noupdate"},description="Don't query the variant already having the tag / do not update the existing annotation")
	private boolean dontUpdateIfInfoIsPresent = false;
	@Parameter(names={"--stopOnError"},description="Stop on network error.")
	private boolean stopOnNetworkError = false;
	@Parameter(names={"--baseurl"},description="Beacon Base URL API")
	private String baseurl="https://beacon-network.org/api";
	@Parameter(names={"--cert"},description="ignore SSL certification errors")
	private boolean ignoreCertErrors = false;
	@Parameter(names={"--tee"},description="show what's happening in the network")
	private boolean teeInput = false;
	
	/** BerkeleyDB Environment to store results */
	private Environment bdbEnv=null;
	/** BerkeleyDB beacon buffer */
	private Database beaconDatabase=null;
	/** BerkeleyDB transaction */
	private Transaction txn=null;
	
	private static class StoredResponse
		{
		long timeStamp;
		Set<String> foundIn=new HashSet<>();
		}
	private static class StoredResponseBinding extends TupleBinding<StoredResponse>
		{
		@Override
		public StoredResponse entryToObject(TupleInput in) {
			StoredResponse st=new StoredResponse();
			st.timeStamp = in.readLong();
			int n= in.readInt();
			for(int i=0;i< n;++i) st.foundIn.add(in.readString());
			return st;
		}
		@Override
		public void objectToEntry(StoredResponse st, TupleOutput out) {
			out.writeLong(st.timeStamp);
			out.writeInt(st.foundIn.size());
			for(final String sw:st.foundIn)
				out.writeString(sw.toString());
			}
		}

	private InputStream wrapTee(final InputStream in) {
		if(this.teeInput) {
			return new TeeInputStream(in, stderr(),false);
			}
		else
			{
			return in;
			}
		}
	
	@Override
	protected int doVcfToVcf(String inputName,final VCFIterator iter,final VariantContextWriter out) {
		CloseableHttpClient httpClient=null;
		InputStream contentInputStream = null;
		
		try 
			{


		   final org.apache.http.impl.client.HttpClientBuilder hb=
				   HttpClients.
				   custom().
				   useSystemProperties()
				   ;
			
			if (this.ignoreCertErrors) {
				// http://stackoverflow.com/questions/24720013/apache-http-client-ssl-certificate-error
				System.setProperty("jsse.enableSNIExtension", "false");
				final SSLContext sslContext =

				org.apache.http.conn.ssl.SSLContexts.custom()
						.loadTrustMaterial(null, new org.apache.http.conn.ssl.TrustStrategy() {
							@Override
							public boolean isTrusted(final X509Certificate[] chain, final String authType)
									throws CertificateException {
								return true;
							}
						}).
						useTLS().
						build();

			final org.apache.http.conn.ssl.SSLConnectionSocketFactory connectionFactory = new org.apache.http.conn.ssl.SSLConnectionSocketFactory(
					sslContext,
					new org.apache.http.conn.ssl.AllowAllHostnameVerifier()
					);

			hb.setSSLSocketFactory(connectionFactory);

			}
			//proxy
			{
	
			String proxyHost = System.getProperty("http.proxyHost","");
			String proxyPortNum = System.getProperty("http.proxyPort","");
			if(!StringUtil.isBlank(proxyHost) && !StringUtil.isBlank(proxyPortNum)) {
				 hb.setProxy(new HttpHost(proxyHost, Integer.parseInt(proxyPortNum),"http"));
				}
			
			proxyHost = System.getProperty("https.proxyHost","");
			proxyPortNum = System.getProperty("https.proxyPort","");
			if(!StringUtil.isBlank(proxyHost) && !StringUtil.isBlank(proxyPortNum)) {
				hb.setProxy(new HttpHost(proxyHost, Integer.parseInt(proxyPortNum),"https"));
				}
			
			
			}
			
			httpClient = hb.build(); 	
			HttpGet httpGetRequest = null;

			final Set<String> available_chromosomes = new HashSet<>();
			try {
				
				httpGetRequest = new HttpGet(baseurl+"/chromosomes");
				httpGetRequest.setHeader("Accept", ContentType.APPLICATION_JSON.getMimeType());
				contentInputStream = httpClient.execute(httpGetRequest).getEntity().getContent();
				final JsonParser jsonparser = new JsonParser();
				final JsonElement root = jsonparser.parse(new InputStreamReader(wrapTee(contentInputStream)));
				Iterator<JsonElement> jsr = root.getAsJsonArray().iterator();
				while (jsr.hasNext()) {
					final String ctg = jsr.next().getAsString();
					available_chromosomes.add(ctg);
				}
				LOG.debug(available_chromosomes);
			} catch (final Exception err) {
				LOG.error(err);
				return -1;
			} finally {
				CloserUtil.close(contentInputStream);
			}

			final Set<String> available_alleles = new HashSet<>();

			try {
				httpGetRequest = new HttpGet(baseurl+"/alleles");
				httpGetRequest.setHeader("Accept", ContentType.APPLICATION_JSON.getMimeType());
				contentInputStream = httpClient.execute(httpGetRequest).getEntity().getContent();

				JsonParser jsonparser = new JsonParser();
				final JsonElement root = jsonparser.parse(new InputStreamReader(contentInputStream));
				Iterator<JsonElement> jsr = root.getAsJsonArray().iterator();
				while (jsr.hasNext()) {
					final String allele = jsr.next().getAsString();
					available_alleles.add(allele);
				}
				LOG.debug(available_alleles);
			} catch (final Exception err) {
				LOG.error(err);
				return -1;
			} finally {
				CloserUtil.close(contentInputStream);
			}

			final StoredResponseBinding storedResponseBinding = new StoredResponseBinding();
			final VCFHeader header = new VCFHeader(iter.getHeader());

			final VCFInfoHeaderLine infoHeaderLine = new VCFInfoHeaderLine(this.infoTag, VCFHeaderLineCount.UNBOUNDED,
					VCFHeaderLineType.String, "Tag inserted with " + getProgramName());
			header.addMetaDataLine(infoHeaderLine);
			DatabaseEntry key = new DatabaseEntry();
			DatabaseEntry data = new DatabaseEntry();
			out.writeHeader(header);
			while (iter.hasNext()) {
				final VariantContext ctx = iter.next();
				if (!ctx.isVariant() || ctx.getReference().isSymbolic()) {
					out.add(ctx);
					continue;
				}

				if (ctx.hasAttribute(infoHeaderLine.getID()) && this.dontUpdateIfInfoIsPresent) {
					out.add(ctx);
					continue;
				}

				String beaconContig = ctx.getContig();
				if (!available_chromosomes.contains(beaconContig)) {
					if (beaconContig.startsWith("chr")) {
						beaconContig = beaconContig.substring(3);
					}
					if (!available_chromosomes.contains(beaconContig)) {
						out.add(ctx);
						continue;
					}

				}

				final List<Allele> altAlleles = ctx.getAlternateAlleles();
				if (altAlleles.isEmpty()) {
					out.add(ctx);
					continue;
				}
				final Set<String> newInfo = new HashSet<>();
				for (final Allele alt : altAlleles) {
					if (alt.isSymbolic() || alt.isNoCall())
						continue;
					final StringBuilder buildUrl = new StringBuilder();
					buildUrl.append("chrom=");
					buildUrl.append(StringUtils.escapeHttp(beaconContig));
					buildUrl.append("&pos=");
					/*
					 * "Coordinate within a chromosome. Position is a number and is 0-based"
					 * .
					 */
					buildUrl.append(ctx.getStart() - 1);
					buildUrl.append("&allele=");

					final String allele;

					if (ctx.getReference().length() > alt.length()) {
						allele = "D";// del
					} else if (ctx.getReference().length() > alt.length()) {
						allele = "I";// ins
					} else {
						allele = alt.getDisplayString();
					}
					if (!available_alleles.contains(allele))
						continue;
					buildUrl.append(allele);
					buildUrl.append("&ref=");
					buildUrl.append(StringUtils.escapeHttp(this.genomeBuild));

					final String queryUrl = buildUrl.toString();

					boolean foundInBdb = false;
					Set<String> foundIn = null;
					if (this.beaconDatabase != null) {
						StringBinding.stringToEntry(queryUrl, key);
						if (this.beaconDatabase.get(this.txn, key, data, LockMode.DEFAULT) == OperationStatus.SUCCESS) {
							StoredResponse response = storedResponseBinding.entryToObject(data);
							if (response.timeStamp < 0) // TODO check how old is
														// that data
							{
								response = null;
								this.beaconDatabase.delete(this.txn, key);
							}
							if (response != null) {
								foundInBdb = true;
								foundIn = response.foundIn;
							}
						}
					}

					if (foundIn == null) {
						foundIn = new HashSet<>();

						
						try {
							httpGetRequest = new HttpGet(baseurl+"/responses?" + queryUrl);
							httpGetRequest.setHeader("Accept", ContentType.APPLICATION_JSON.getMimeType());
							
							LOG.debug(httpGetRequest.getURI());
							
							
							contentInputStream = httpClient.execute(httpGetRequest).getEntity().getContent();

							final JsonParser jsonparser = new JsonParser();
							final JsonElement root = jsonparser.parse(new InputStreamReader(contentInputStream));
							final Iterator<JsonElement> jsr = root.getAsJsonArray().iterator();
							while (jsr.hasNext()) {
								final JsonObject b = jsr.next().getAsJsonObject();
								if (!(b.has("beacon") && b.has("response")))
									continue;
								final String beacon_id = b.get("beacon").getAsJsonObject().get("id").getAsString();
								final JsonElement response_prim = b.get("response");
								if (response_prim.isJsonPrimitive() && response_prim.getAsBoolean()) {
									foundIn.add(beacon_id);
								}
							}

						} catch (final Exception err) {
							if (stopOnNetworkError) {
								LOG.error(err);
								throw new RuntimeIOException(err);
								}
							else
								{
								LOG.warn(err.getMessage());
								}
						}
						finally {
							CloserUtil.close(contentInputStream);
							}
						}
					

					if (this.beaconDatabase != null && !foundInBdb) {
						StoredResponse response = new StoredResponse();
						response.timeStamp = System.currentTimeMillis();
						response.foundIn = foundIn;
					}
					// 17&pos=41244981&=G&ref=GRCh37")
					newInfo.addAll(
							foundIn.stream().map(S -> alt.getDisplayString() + "|" + S).collect(Collectors.toSet()));
				}
				if (newInfo.isEmpty()) {
					out.add(ctx);
					continue;
				}

				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				vcb.attribute(infoHeaderLine.getID(), new ArrayList<String>(newInfo));
				out.add(vcb.make());
			}
			return 0;
		}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(httpClient);
			}
		}
	
	@Override
		public int doWork(final List<String> args) {
			try
				{ 
				
				
				if(this.bdbDir!=null) {
					LOG.info("open BDB "+this.bdbDir);
					IOUtil.assertDirectoryIsWritable(this.bdbDir);
					final EnvironmentConfig envCfg=new EnvironmentConfig();
					envCfg.setAllowCreate(true);
					envCfg.setReadOnly(false);
					this.bdbEnv = new Environment(this.bdbDir, envCfg);
					
					final DatabaseConfig cfg=new DatabaseConfig();
					cfg.setAllowCreate(true);
					cfg.setReadOnly(false);
					this.beaconDatabase = this.bdbEnv.openDatabase(this.txn,"ga4ghBeaconBuffer",cfg);
					}
				return doVcfToVcf(args, outputFile);
				}
			catch(final Exception err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(beaconDatabase);
				CloserUtil.close(bdbEnv);
				}
			
			}
	
public static void main(final String[] args) {
	new VcfAnnotWithBeacon().instanceMainWithExit(args);
	}
}
