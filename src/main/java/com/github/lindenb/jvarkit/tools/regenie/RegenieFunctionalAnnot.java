package com.github.lindenb.jvarkit.tools.regenie;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
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
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.ucsc.UcscTranscript;
import com.github.lindenb.jvarkit.ucsc.UcscTranscriptReader;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

/**
BEGIN_DOC

## Example

snpeff.cadd.in.vcf.gz |\
	java -jar dist/jvarkit.jar regeniefunctionalannot -A anot.tsv |\
	java -jar dist/jvarkit.jar regeniemakeannot -o OUT


END_DOC
*/
@Program(name="regeniefunctionalannot",
	description="Create annotation files for regenie using snpEff annotations",
	keywords={"vcf","regenie","burden"},
	creationDate="20250311",
	modificationDate="20250320",
	jvarkit_amalgamion = true,
	generate_doc = true
	)
public class RegenieFunctionalAnnot extends AbstractRegenieAnnot {
	static final String FIRST_INTRON="first_intron";
	private static final Logger LOG = Logger.of(RegenieFunctionalAnnot.class);
	@Parameter(names={"-A","--annotations"},description="seq_ontology <-> score file. TSV file. no header. at least 2 columns prediction_name/score",required = true)
	private Path masksFile = null;
	@Parameter(names={"--kg","--known"},description="known gene data for first intron/intergenic. " + UcscTranscriptReader.OPT_DESC)
	private Path knownGene = null;

	
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
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		
		if(knownGene!=null && !prediction2score.containsKey(FIRST_INTRON)) {
			LOG.warning("No "+FIRST_INTRON+" defined in "+masksFile);
		}
		
		if(knownGene!=null && prediction2score.containsKey(FIRST_INTRON)) {
			this.firstIntronToTranscript = new IntervalTreeMap<>();
			try(UcscTranscriptReader rd=new UcscTranscriptReader(this.knownGene)) {
				try(CloseableIterator<UcscTranscript> iter= rd.iterator()) {
					while(iter.hasNext()) {
						final UcscTranscript kg = iter.next();
						if(kg.getIntronCount()<=0) continue;
						final String ctg = fixContig(kg.getContig());
						if(StringUtils.isBlank(ctg)) continue;
						final	UcscTranscript.Intron intron;
						if(kg.isPositiveStrand()) {
							intron = kg.getIntron(0);
							}
						else
							{
							intron = kg.getIntron(kg.getIntronCount()-1);
							}
						final Interval r = new Interval(ctg,intron.getStart(),intron.getEnd());
						Set<String> transcriptids = firstIntronToTranscript.get(r);
						if(transcriptids==null) {
							transcriptids=new HashSet<>();
							firstIntronToTranscript.put(r, transcriptids);
							}
						transcriptids.add(kg.getTranscriptId());
						}
					}
				catch(IOException err) {
					throw new RuntimeIOException(err);
					}
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
						if (score == null) throw new IOException("undefined prediction key "+pred_key+ " in file "+masksFile+" available are:"+String.join(",", this.prediction2score.keySet()));
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
					v.score = OptionalDouble.of( best_score);
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
