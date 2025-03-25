package com.github.lindenb.jvarkit.tools.regenie;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;

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
	private static final Logger LOG = Logger.build(RegenieFunctionalAnnot.class).make();
	@Parameter(names={"-A","--annotations"},description="seq_ontology <-> score file. TSV file. no header. at least 2 columns prediction_name/score",required = true)
	private Path masksFile = null;

	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	private AnnPredictionParser annParser;
	private final Map<String,Double> prediction2score= new HashMap<>();
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
					v.contig = fixContig(ctx.getContig());
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
		}

	
	

	public static void main(String[] args) {
		new RegenieFunctionalAnnot().instanceMainWithExit(args);
	}

}
