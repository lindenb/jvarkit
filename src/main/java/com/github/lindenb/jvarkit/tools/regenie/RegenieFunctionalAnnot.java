/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.regenie;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.gtf.GTFCodec;
import com.github.lindenb.jvarkit.gtf.GTFLine;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

/**
BEGIN_DOC

 The aim of this  class is to produce a file for regenie containing functional annotations with the following header:
 
<pre>"CONTIG","POS","ID","GENE","ANNOTATION","SCORE","CADD","FREQ","SINGLETON"</pre>
 


## Example


END_DOC
*/
@Program(name="regeniefunctionalannot",
	description="Create annotation files for regenie using snpEff annotations",
	keywords={"vcf","regenie","burden"},
	creationDate="20250311",
	modificationDate="20260315",
	jvarkit_amalgamion = true,
	generate_doc = true
	)
public class RegenieFunctionalAnnot extends AbstractRegenieAnnot {
	static final String FIRST_INTRON="first_intron";
	private static final Logger LOG = Logger.of(RegenieFunctionalAnnot.class);
	@Parameter(names={"-A","--annotations"},description="seq_ontology <-> score file. TSV file. no header. at least 2 columns prediction_name/score",required = true)
	private Path masksFile = null;
	@Parameter(names={"--gtf"},description="GTF file used to get first intron/intergenic. ")
	private Path gtfPath = null;

	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	private AnnPredictionParser annParser;
	private final Map<String,Double> prediction2score= new HashMap<>();
	private IntervalTreeMap<Set<String>> firstIntronToTranscript=null;
	@Override
	protected VCFHeader initVcfHeader(VCFHeader h) {
		h =  super.initVcfHeader(h);
		this.annParser = new AnnPredictionParserFactory(h).get();
		if(!this.annParser.isValid()) throw new IllegalArgumentException("cannot create ANN parser");
		
		try(BufferedReader br= IOUtils.openPathForBufferedReading(masksFile)) {
			final Pattern spaces_regex = Pattern.compile("[ \t]");
				for(;;) {
					String line = br.readLine();
					if(line==null) break;
					if(line.startsWith("#")) continue;
					
					final String[] tokens = spaces_regex.split(line);
					if(tokens.length<2) throw new JvarkitException.TokenErrors(2, tokens);
					final String pred = tokens[0];
					if(StringUtils.isBlank(pred)) {
						throw new IllegalArgumentException("empty prediction in "+line+" in "+masksFile);
						}
					if(prediction2score.containsKey(tokens[0])) {
						throw new IllegalArgumentException("duplicate prediction "+tokens[0]+" in "+masksFile);
						}
					prediction2score.put(pred, Double.parseDouble(tokens[1].trim()));
		            }
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		
		if(gtfPath!=null && !prediction2score.containsKey(FIRST_INTRON)) {
			LOG.warning("No "+FIRST_INTRON+" defined in "+masksFile);
			}
		
		if(gtfPath!=null && prediction2score.containsKey(FIRST_INTRON)) {
			this.firstIntronToTranscript = new IntervalTreeMap<>();
			try(BufferedReader rd= IOUtils.openPathForBufferedReading(this.gtfPath)) {
				final GTFCodec gtfCodec = new GTFCodec();
				final Map<String,List<Interval>> tr2exons = new HashMap<>();
				final SAMSequenceDictionary dict = h.getSequenceDictionary();
				if(dict!=null) {
					gtfCodec.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
					}
				for(;;) {
					final String line = rd.readLine();
					if(line==null) break;
					if(line.startsWith("#") || StringUtils.isBlank(line)) continue; 
					final GTFLine rec= gtfCodec.decode(line);
					if(rec==null) continue;
					if(!rec.getType().equals("exon")) continue;
					if(!rec.hasAttribute("transcript_id")) continue;
					final String tr = removeVersionFromEnst(rec.getAttribute("transcript_id"));
					List<Interval> exons = tr2exons.get(tr);
					if(exons==null) {
						exons = new ArrayList<>();
						tr2exons.put(tr, exons);
						}
					exons.add(new Interval(rec.getContig(),rec.getStart(),rec.getEnd(),rec.isNegativeStrand(),tr));
					}
				for(String tr:tr2exons.keySet()) {
					List<Interval> L = tr2exons.get(tr);
					if(L.size()<2) continue;//no intron
					Collections.sort(L,(A,B)->Integer.compare(A.getStart(), B.getStart()));
					final Interval first = L.get(0);
					final Interval intron;
					if(first.isPositiveStrand()) {
						intron = new Interval(
								first.getContig(),
								L.get(0).getEnd()+1,
								L.get(1).getStart()-1,
								false,
								first.getName()
								);
						}
					else {
						intron = new Interval(
								first.getContig(),
								L.get(L.size()-2).getEnd()+1,
								L.get(L.size()-1).getStart()-1,
								false,
								first.getName()
								);
						}
					Set<String> transcriptids = firstIntronToTranscript.get(intron);
					if(transcriptids==null) {
						transcriptids=new HashSet<>();
						firstIntronToTranscript.put(intron, transcriptids);
						}
					transcriptids.add(first.getName());
					}
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		
		return h;
		}
	
	private String removeVersionFromEnst(final String id) {
		if(id.startsWith("ENST") ) {
			int dot = id.lastIndexOf(".");
			if(dot>0) return id.substring(0,dot);
		}
		return id;
		}
	
	@Override
	protected void dump(final PrintWriter w,final VariantContext ctx) throws Exception {
		final String contig = fixContig(ctx.getContig());
		final String altstr = ctx.getAlternateAllele(0).getDisplayString();
		final List<AnnPredictionParser.AnnPrediction> predictions = this.annParser.getPredictions(ctx);
		
		
		// two loop, one for Transcript, the other for gene
		for(int side=0;side< 2;++side) {
			final Function<AnnPredictionParser.AnnPrediction, String> extract_gene = side==0?
					PRED->PRED.getGeneName():
					PRED->removeVersionFromEnst(PRED.getFeatureId())
					;
			
			final Set<String> gene_names=predictions.stream().
					map(extract_gene).
					filter(S->!(S.isEmpty()|| S.equals("."))).
					collect(Collectors.toSet());
					
			
			for(String gene_name:gene_names) {
				Double best_score = null;
				String best_pred = null;
				for(AnnPredictionParser.AnnPrediction pred:predictions) {
					if(!gene_name.equals(extract_gene.apply(pred))) continue;
					if(!pred.getAllele().equalsIgnoreCase(altstr)) continue;
					for(String pred_key : pred.getSOTermsStrings()) {
						if(pred_key.equals("intergenic_region")) continue;
						final Double score = this.prediction2score.getOrDefault(pred_key, null);
						if (score == null) {
							throw new IOException("undefined prediction key "+pred_key+ " in file "+masksFile+" available are:"+String.join(",", this.prediction2score.keySet()));
							}
						if (best_score == null || best_score.compareTo(score) < 0) {
							best_score =score;
							best_pred = pred_key;
							}
						}
					}
				
				if (best_score != null) {
					final Variation v = new Variation();
					v.contig = contig;
					v.pos = ctx.getStart();
					v.id = makeID(ctx);
					v.gene = gene_name;
					v.prediction = best_pred;
					if(super.isIgnoringMaskScore()) {
							v.score = OptionalDouble.of(1.0);
							}
						else
							{
                            v.score = OptionalDouble.of( best_score);
                            }
					v.cadd = getCaddScore(ctx);
					v.is_singleton = isSingleton(ctx);
					v.frequency = getFrequency(ctx);
					print(w,v);
					}
				}
			}
		
		if(this.firstIntronToTranscript!=null) {
			final Double score = this.prediction2score.getOrDefault(FIRST_INTRON, null);
			for(String transcriptId : this.firstIntronToTranscript.
					getOverlapping(new Interval(contig,ctx.getStart(),ctx.getEnd())).stream().
					flatMap(SET->SET.stream()).
					collect(Collectors.toSet())) {
					final Variation v = new Variation();
					v.contig = contig;
					v.pos = ctx.getStart();
					v.id = makeID(ctx);
					v.gene = transcriptId;
					v.prediction = FIRST_INTRON;
					v.score = OptionalDouble.of( score);
					v.cadd = getCaddScore(ctx);
					v.is_singleton = isSingleton(ctx);
					v.frequency = getFrequency(ctx);
					print(w,v);
					}
				}
		}

	
	

	public static void main(String[] args) {
		new RegenieFunctionalAnnot().instanceMainWithExit(args);
	}

}
