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
package com.github.lindenb.jvarkit.tools.gtf;

import java.io.BufferedReader;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.nio.file.Path;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.tribble.TribbleException;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3Constants;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.gff.Gff3FeatureImpl;
import htsjdk.tribble.gff.Gff3Writer;

/**

BEGIN_DOC

```
```

END_DOC
*/
@Program(
		name="gtf2gff",
		description="Convert GTF to gff",
		creationDate="20220703",
		modificationDate="20250328",
		keywords= {"gtf","gff","gff3"}
		)
public class GtfToGff
	extends Launcher {
	private static final Logger LOG = Logger.build(GtfToGff.class).make();
	
	@Parameter(names={"-o","--output"},description= OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--biotype"},description= "Default biotype (empty:no default)")
	private String default_biotype=null;
	@Parameter(names={"--make-gene"},description= "When 'gene' is missing, create a gene for each 'transcript/mrna'")
	private boolean make_gene = false;
	@Parameter(names={"--escape-UTF8"},description= "escape UTF-8")
	private boolean escape_UTF8 = false;


	private String escape(final String s) {
		if(!escape_UTF8) return s;
		 try {
	            return URLEncoder.encode(s, "UTF-8").replace("+", " ");
	        } catch (final UnsupportedEncodingException ex) {
	            throw new TribbleException("Encoding failure", ex);
	        }
	}

	@Override
	public int doWork(final List<String> args)  {
		try {
			final Predicate<String> isTranscript =  S-> 
				S.equalsIgnoreCase("transcript") ||
				S.equalsIgnoreCase("mrna")
				;
			final Predicate<String> isChildOfTranscript =  S-> 
				S.equalsIgnoreCase("exon") ||
				S.equalsIgnoreCase("CDS") ||
				S.equalsIgnoreCase("three_prime_UTR") || 
				S.equalsIgnoreCase("five_prime_UTR") || 
				S.equalsIgnoreCase("3_prime_UTR") || 
				S.equalsIgnoreCase("5_prime_UTR") || 
				S.equalsIgnoreCase("UTR") || 
				S.equalsIgnoreCase("start_codon") ||
				S.equalsIgnoreCase("stop_codon")
				;
			final String input = oneFileOrNull(args);
				try(final BufferedReader br = super.openBufferedReader(input)) {
				try(Gff3Writer out = this.outputFile==null?
							new Gff3Writer(stdout()) {
								protected String escapeString(String s) { return GtfToGff.this.escape(s);}
								}:
							new Gff3Writer(this.outputFile) {
								protected String escapeString(String s) { return GtfToGff.this.escape(s);}
								}
							) {
					String line;
					final GTFCodec codec = new GTFCodec();
					while((line=br.readLine())!=null) {
						if(StringUtils.isBlank(line)) continue;
						if(line.startsWith("#")) {
							out.addComment(line.substring(1));
							continue;
							}
						final String[] tokens = CharSplitter.TAB.split(line);
						for(int side=0;side< (make_gene?2:1);++side) {
							if(side==2) {
								if(isTranscript.test(tokens[2])) {
									tokens[2]="gene";
									}
								else {
									break;
									}
								}
							
							final GTFLine rec = codec.decode(tokens);
							if(rec==null) continue;
							final  Map<String, String> hash = rec.getAttributes();
							final  Map<String, List<String>> attributes = new LinkedHashMap<>(hash.size());
							
							for(String key: hash.keySet()) {
								final List<String> value =Collections.singletonList(hash.get(key));
								attributes.put(key, value);
								}
							
							if(rec.getType().equals("gene")  && rec.hasAttribute("gene_id")) {
								attributes.put(Gff3Constants.ID_ATTRIBUTE_KEY, Collections.singletonList("gene:"+rec.getAttribute("gene_id")));
								if(!attributes.containsKey(Gff3Constants.NAME_ATTRIBUTE_KEY)) {
									attributes.put(Gff3Constants.NAME_ATTRIBUTE_KEY, Collections.singletonList(rec.getAttribute("gene_id")));
									if(rec.hasAttribute("gene_name")) {
										attributes.put(Gff3Constants.NAME_ATTRIBUTE_KEY, Collections.singletonList(rec.getAttribute("gene_name")));
										}
									}
								attributes.remove("transcript_id");
								attributes.remove("transcript_name");
								}
							else if(isTranscript.test(rec.getType())) {
								if(!StringUtils.isBlank(this.default_biotype)) {
									if(!(attributes.containsKey("biotype")||  attributes.containsKey("transcript_biotype") )) {
										attributes.put("transcript_biotype", Collections.singletonList(this.default_biotype));
										attributes.put("biotype", Collections.singletonList(this.default_biotype));
										}
									}
								if(rec.hasAttribute("transcript_id")) {
									attributes.put(Gff3Constants.ID_ATTRIBUTE_KEY, Collections.singletonList("transcript:"+rec.getAttribute("transcript_id")));
									}
								if(rec.hasAttribute("gene_id")) {
									attributes.put(Gff3Constants.PARENT_ATTRIBUTE_KEY, Collections.singletonList("gene:"+rec.getAttribute("gene_id")));
									}
								if(!attributes.containsKey(Gff3Constants.NAME_ATTRIBUTE_KEY)) {
									attributes.put(Gff3Constants.NAME_ATTRIBUTE_KEY, Collections.singletonList(rec.getAttribute("transcript_id")));
									if(rec.hasAttribute("transcript_name")) {
										attributes.put(Gff3Constants.NAME_ATTRIBUTE_KEY, Collections.singletonList(rec.getAttribute("transcript_name")));
										}
									}
								}
							else if(isChildOfTranscript.test(rec.getTranscriptId()) || rec.hasAttribute("transcript_id")) {
								attributes.put(Gff3Constants.ID_ATTRIBUTE_KEY,Collections.singletonList(StringUtils.md5(line).substring(0, 9)));
								if(rec.hasAttribute("transcript_id")) {
									attributes.put(Gff3Constants.PARENT_ATTRIBUTE_KEY, Collections.singletonList("transcript:"+rec.getAttribute("transcript_id")));
									}
								}
							
							final Gff3Feature feat = new Gff3FeatureImpl(
									rec.getContig(),
									rec.getSource(),
									rec.getType(),
									rec.getStart(),
									rec.getEnd(),
									(rec.getScore()==null?-1:rec.getScore()),
									Strand.decode(rec.getStrand()),
									rec.getPhase(),
									attributes
									);
							
							
							out.addFeature(feat);
							}
						}
					}
				}
			
		
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	
	public static void main(final String[] args)
		{
		new GtfToGff().instanceMainWithExit(args);
		}
	}
