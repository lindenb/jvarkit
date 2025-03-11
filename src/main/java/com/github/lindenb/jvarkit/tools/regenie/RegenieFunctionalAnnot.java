package com.github.lindenb.jvarkit.tools.regenie;

import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;

import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

/** 
BEGIN_DOC



END_DOC
*/
@Program(name="regeniefunctionalannot",
	description="Create annotation files for regenie using snpEff annotations",
	keywords={"vcf","regenie","burden"},
	creationDate="20250311",
	modificationDate="20250311"
	)
public class RegenieFunctionalAnnot extends AbstractRegenieAnnot {
	private static final Logger LOG = Logger.build(RegenieFunctionalAnnot.class).make();
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	private AnnPredictionParser annParser;
	
	
	private void initScores() {
		makeScore("3_prime_UTR_variant", 0.1,"UTR", "UTR3");
		makeScore("5_prime_UTR_premature_start_codon_gain_variant", 0.2,"UTR", "UTR5");
		makeScore("5_prime_UTR_truncation", 0.5,"UTR", "UTR5");
		makeScore("3_prime_UTR_truncation", 0.5,"UTR", "UTR3");
		makeScore("5_prime_UTR_variant", 0.2,"UTR", "UTR5");
		makeScore("bidirectional_gene_fusion", 1.0);
		makeScore("conservative_inframe_deletion", 0.1,"protein_altering");
		makeScore("conservative_inframe_insertion", 0.3,"protein_altering");;
		makeScore("disruptive_inframe_deletion", 0.2,"protein_altering");
		makeScore("disruptive_inframe_insertion", 0.2,"protein_altering");
		makeScore("downstream_gene_variant", 0.1,"downstream", "updownstream");
		makeScore("exon_loss_variant", 1.0,"protein_altering");
		makeScore("frameshift_variant", 0.4);
		makeScore("exon_loss_variant", 1.0,"protein_altering");
		makeScore("gene_fusion", 0.9,"protein_altering");
		makeScore("intergenic_region", 0.001);
		makeScore("intragenic_variant", 0.01);
		makeScore("initiator_codon_variant",0.5,"protein_altering");
		makeScore("intron_variant", 0.05,"intronic");
		makeScore("missense_variant", 0.9,"protein_altering");
		makeScore("non_coding_transcript_exon_variant", 0.1,"non_coding");
		makeScore("non_coding_transcript_variant", 0.1,"non_coding");
		makeScore("splice_acceptor_variant", 0.5,"protein_altering", "splice");
		makeScore("splice_donor_variant", 0.5,"protein_altering", "splice");
		makeScore("splice_region_variant", 0.5,"protein_altering", "splice");
		makeScore("start_retained_variant",0.1,"synonymous");
		makeScore("start_lost", 0.6,"protein_altering");
		makeScore("stop_gained", 0.9,"protein_altering");
		makeScore("stop_lost", 0.6,"protein_altering");
		makeScore("stop_retained_variant", 0.2,"synonymous");
		makeScore("synonymous_variant", 0.1,"synonymous");
		makeScore("upstream_gene_variant", 0.1,"upstream", "updownstream");
		}
	
	
	@Override
	protected VCFHeader initVcfHeader(VCFHeader h) {
		initScores();
		h =  super.initVcfHeader(h);
		this.annParser = new AnnPredictionParserFactory(h).get();
		if(!this.annParser.isValid()) throw new IllegalArgumentException("cannot create ANN parser");
		return h;
		}
		
	@Override
	protected void dump(final SortingCollection<Variation> sorter,final VariantContext ctx) throws Exception {
		final String altstr = ctx.getAlternateAllele(0).getDisplayString();
		final List<AnnPredictionParser.AnnPrediction> predictions = this.annParser.getPredictions(ctx);
			
		
		
		for(int side=0;side< 2;++side) {
			final Function<AnnPredictionParser.AnnPrediction, String> extract_gene = side==0?
					PRED->PRED.getGeneName():
					PRED->PRED.getFeatureId()
					;
			
			final Set<String> gene_names=predictions.stream().
					map(extract_gene).
					filter(S->!(S.isEmpty()|| S.equals("."))).
					collect(Collectors.toSet());
					
			
			for(String gene_name:gene_names) {
				Double best_score = null;
				Prediction best_pred = null;
				for(AnnPredictionParser.AnnPrediction pred:predictions) {
					if(!gene_name.equals(extract_gene.apply(pred))) continue;
					if(!pred.getAllele().equalsIgnoreCase(altstr)) continue;
					for(String pred_key : pred.getSOTermsStrings()) {
						if(pred_key.equals("intergenic_region")) continue;
						final Prediction p = super.name2prediction.getOrDefault(pred_key, null);
						if (p == null) throw new IOException("undefined prediction key "+pred_key+ " in "+String.join(",",super.name2prediction.keySet()));
						if (best_score == null || best_score.compareTo(p.score) < 0) {
							best_score = p.score;
							best_pred = p;
							}
						}
					}
				
				if (best_score != null) {
					final Variation v = new Variation();
					v.contig = fixContig(ctx.getContig());
					v.pos = ctx.getStart();
					v.id = ctx.getID();
					v.gene = gene_name;
					v.prediction = best_pred.name;
					v.score = best_score;
					v.cadd = getCaddScore(ctx);
					sorter.add(v);
					}
				}
			}		
		}

	
	

	public static void main(String[] args) {
		new RegenieFunctionalAnnot().instanceMainWithExit(args);
	}

}
