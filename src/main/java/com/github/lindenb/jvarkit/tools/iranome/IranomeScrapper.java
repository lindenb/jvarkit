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
package com.github.lindenb.jvarkit.tools.iranome;


import java.io.InputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import org.apache.http.StatusLine;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

## Example

```
$ java -jar dist/iranomescrapper.jar src/test/resources/test_vcf01.vcf
[WARN][IranomeScrapper]Too many requests for http://www.iranome.ir/variant/1-912049-T-C max_sleep_millisec=1
[WARN][IranomeScrapper]Too many requests for http://www.iranome.ir/variant/1-913889-G-A max_sleep_millisec=1001
CHROM	POS	REF	ALT	AC	AN	AF	AC_ARAB	AC_AZERI	AC_BALOCH	AC_KURD	AC_LUR	AC_PERSIAN	AC_PERSIAN_GULF_ISLANDER	AC_TURKMEN	AN_ARAB	AN_AZERI	AN_BALOCH	AN_KURD	AN_LUR	AN_PERSIAN	AN_PERSIAN_GULF_ISLANDER	AN_TURKMEN	date
1	909238	G	C	1241	1586	0.782472	160	149	171	151	156	158	157	139	200	198	200	200	200	200	200	188	20210629_103622
1	912049	T	C	856	1558	0.549422	104	106	120	85	105	117	116	103	190	196	194	194	196	200	194	194	20210629_103623
1	914333	C	G	961	1584	0.606692	119	117	132	107	121	127	129	109	200	200	200	200	200	200	200	184	20210629_103626
1	914852	G	C	939	1596	0.588346	117	116	131	99	112	122	127	115	200	200	200	200	200	200	200	196	20210629_103628
1	914940	T	C	937	1598	0.586358	116	116	131	99	112	122	127	114	200	200	200	200	200	200	200	198	20210629_103628
1	935222	C	A	989	1554	0.636422	123	121	139	108	127	130	139	102	200	198	200	200	200	200	200	156	20210629_103638
1	949608	G	A	460	1600	0.2875	55	59	60	62	50	59	56	59	200	200	200	200	200	200	200	200	20210629_103647
1	984302	T	C	934	1554	0.60103	116	120	122	115	125	124	125	87	200	200	200	198	200	200	200	156	20210629_103719
1	985266	C	T	934	1600	0.58375	115	118	119	108	126	117	124	107	200	200	200	200	200	200	200	200	20210629_103720
1	985446	G	T	647	1330	0.486466	79	91	97	70	81	87	99	43	166	182	186	176	168	158	174	120	20210629_103721
1	985450	G	A	257	1352	0.190089	33	36	42	24	30	24	38	30	166	176	180	176	172	172	176	134	20210629_103723
1	1007432	G	A	908	1598	0.56821	114	111	109	115	113	115	123	108	200	200	200	200	200	200	200	198	20210629_103735
1	1018144	T	C	907	1596	0.568296	115	111	108	115	113	115	127	103	200	200	200	200	200	200	200	196	20210629_103746
1	1019180	T	C	947	1600	0.591875	116	114	113	118	118	119	131	118	200	200	200	200	200	200	200	200	20210629_103750
```

END_DOC
*/

@Program(name="iranomescrapper",
description="Iranome scrapper",
keywords={"iranome","vcf"},
creationDate="20210728",
modificationDate="20210729"
)
public class IranomeScrapper extends Launcher {
	private static final Logger LOG = Logger.build(IranomeScrapper.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-T","--tabix"},description="Existing tabix-indexed database. Don't lookup the data if the variant is already in this database.")
	private String tabixDatabase = null;
	@Parameter(names={"-failed","--failed"},description="Save the variants where no data in iranome data was found in this file.")
	private Path failedVcfPath = null;

	
	private CloseableHttpClient httpClient = null;
	private final List<String> fields = Arrays.asList("chrom","pos","ref","alt","allele_count","allele_num","allele_freq");
	private final List<String> populations = Arrays.asList("Arab","Azeri","Baloch","Kurd","Lur","Persian","Persian Gulf Islander","Turkmen");
	private long lastCall = System.currentTimeMillis();
	private long max_sleep_millisec = 1;

	
	private static String toIranomeContig(final String c) {
		if(c.startsWith("chr")) return c.substring(3);
		return c;
		}
	
	private Map<String,Object> getIranomeVariant(final VariantContext ctx,final Allele alt) {
		final String iranomeCtg = toIranomeContig(ctx.getContig());
		if(StringUtils.isBlank(iranomeCtg)) return null;
		final Map<String,Object> hash = new HashMap<String, Object>();
		final String urlStr = "http://www.iranome.ir/variant/"+iranomeCtg+"-"+ctx.getStart()+"-"+ctx.getReference().getDisplayString()+"-"+alt.getDisplayString();
		final JsonParser jsonparser = new JsonParser();
		
		try {
			String html = "";
			final HttpGet httpGetRequest = new HttpGet(urlStr);
			for(;;) {
				final long now = System.currentTimeMillis();
				final long waitMillisec = now-this.lastCall;
				this.lastCall = now;
				if(waitMillisec<  this.max_sleep_millisec) {
					TimeUnit.MILLISECONDS.sleep(this.max_sleep_millisec-waitMillisec);
					}
				
				try(CloseableHttpResponse resp = this.httpClient.execute(httpGetRequest)) {
					final StatusLine status = resp.getStatusLine();
					if(status.getStatusCode()!=200) {
						if(status.getStatusCode()==400) return null;
						if(status.getStatusCode()==429) {
							LOG.warning("Too many requests for "+urlStr+" max_sleep_millisec="+max_sleep_millisec);
							this.max_sleep_millisec += 1_000;
							continue;
							}
						LOG.error("cannot get "+urlStr+" "+status.getReasonPhrase()+" ("+status.getStatusCode()+")");
						return null;
						}
					try(InputStream contentInputStream = resp.getEntity().getContent()) {
						html = IOUtil.slurp(contentInputStream);				
						}
					}
				break;
				}
			final String window_variant= "window.variant = ";
			final String end_json = "};";
			final int i1= html.indexOf(window_variant);
			if(i1==-1) return null;
			int i2 = html.indexOf(end_json);
			if(i1==-1) return null;
			final String jsonString = html.substring(i1+window_variant.length(),i2+1).trim();
			
			final JsonObject jobj = jsonparser.parse(jsonString).getAsJsonObject();
			for(final String key:fields) {
				if(!jobj.has(key)) {
					LOG.info("no "+key +" in "+jobj+" "+urlStr);
					return null;
					}
				hash.put(key, jobj.get(key).getAsString());
				}
			if(!jobj.has("pop_acs")|| !jobj.get("pop_acs").isJsonObject()) {
				LOG.error("no pop_acs");
				return null;
				}
			final JsonObject pop_acs = jobj.get("pop_acs").getAsJsonObject();
			for(final String pop:populations) {
				if(!pop_acs.has(pop)) {
					LOG.info("no "+ pop +" in "+pop_acs+" "+urlStr);
					return null;
					}
				hash.put("AC_"+pop, pop_acs.get(pop).getAsInt());
				}
			
			if(!jobj.has("pop_ans") || !jobj.get("pop_ans").isJsonObject()) {
				LOG.error("no pop_acs");
				return null;
				}
			final JsonObject pop_ans = jobj.get("pop_ans").getAsJsonObject();
			for(final String pop:populations) {
				if(!pop_ans.has(pop)) {
					LOG.info("no "+ pop +" in "+pop_ans+" "+urlStr);
					return null;
					}
				hash.put("AN_"+pop, pop_ans.get(pop).getAsInt());
				}
			return hash;
			}
		catch (final Throwable err) {
			LOG.error(err);
			return null;
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		TabixReader tabixReader=null;
		VariantContextWriter failedWriter =null;
		try {
			/** create http client */
			this.httpClient = HttpClients.createSystem();//createDefault();
		
			try(VCFIterator iter = super.openVCFIterator(oneFileOrNull(args))) {
				if(this.tabixDatabase!=null) {
					tabixReader = new TabixReader(this.tabixDatabase);
					}
				
				if(this.failedVcfPath!=null)  {
					failedWriter = VCFUtils.createVariantContextWriterToPath(failedVcfPath);
					failedWriter.writeHeader(iter.getHeader());
					}
				
				try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
					pw.print("CHROM\tPOS\tREF\tALT\tAC\tAN\tAF");
					for(int i=0;i< populations.size();i++) {
						pw.print("\tAC_"+populations.get(i).toUpperCase().replace(' ', '_'));
						}
					for(int i=0;i< populations.size();i++) {
						pw.print("\tAN_"+populations.get(i).toUpperCase().replace(' ', '_'));
						}
					pw.println("\tdate");
					while(iter.hasNext()) {
						final VariantContext ctx = iter.next();
						if(!AcidNucleics.isATGC(ctx.getReference())) continue;
						final String iranomeCtg = toIranomeContig(ctx.getContig());
						for(final Allele alt: ctx.getAlternateAlleles()) {
							if(!AcidNucleics.isATGC(alt)) continue;
							if(tabixReader!=null) {
								boolean found=false;
								final TabixReader.Iterator iter2 = tabixReader.query(iranomeCtg, ctx.getStart(), ctx.getEnd());
								for(;;) {
									final String line2 = iter2.next();
									if(line2==null) break;
									final String[] tokens = CharSplitter.TAB.split(line2);
									if(tokens.length<4) throw new JvarkitException.TokenErrors(4, tokens);
									if(tokens[0].equals(iranomeCtg) &&
										Integer.parseInt(tokens[1])==ctx.getStart() &&
										tokens[2].equals(ctx.getReference().getDisplayString()) &&
										tokens[3].equals(alt.getDisplayString())) {
										LOG.info("already in database "+String.join("-", Arrays.asList(tokens).subList(0,4)));
										found = true;
										break;
									}
								if(found) continue;
								}//end tabix terator
							}//end tabix!=null
						final Map<String, Object> hash = this.getIranomeVariant(ctx,alt);
						if(hash==null) {
							if(failedWriter!=null) {
								failedWriter.add(ctx);
								}
							continue;
							}
						for(int i=0;i< fields.size();i++) {
							pw.print(i>0?"\t":"");
							pw.print(hash.get(fields.get(i)));
							}
						for(int i=0;i< populations.size();i++) {
							pw.print("\t"+hash.get("AC_"+populations.get(i)));
							}
						for(int i=0;i< populations.size();i++) {
							pw.print("\t"+hash.get("AN_"+populations.get(i)));
							}
						pw.print("\t");
						pw.println(StringUtils.now());
						}
					}
				pw.flush();
				}// end print
			}// end iter
			return 0;
		} catch (final Throwable err) {
			LOG.error(err);
			return -1;
		} finally {
			CloserUtil.close(tabixReader);
			CloserUtil.close(failedWriter);
			CloserUtil.close(this.httpClient);
			}
		}
	
	public static void main(String[] args) {
		new IranomeScrapper().instanceMainWithExit(args);
		}
	}
