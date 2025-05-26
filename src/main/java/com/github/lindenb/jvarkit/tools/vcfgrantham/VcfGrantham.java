/**
 * 
 */
package com.github.lindenb.jvarkit.tools.vcfgrantham;

import java.util.Arrays;
import java.util.List;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.bio.AminoAcids;
import com.github.lindenb.jvarkit.bio.AminoAcids.AminoAcid;
import com.github.lindenb.jvarkit.grantham.GranthamScore;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.BcfToolsPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.BcfToolsPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

## Example

```
$ java -jar dist/jeter.jar src/test/resources/rotavirus_rf.ann.vcf.gz | grep GRANT -m3
##INFO=<ID=GRANTHAM_SCORE,Number=1,Type=Integer,Description="Max grantham score for this variant">
RF01	970	.	A	C	48.67	.	AC=2;AN=10;ANN=C|missense_variant|MODERATE|Gene_18_3284|Gene_18_3284|transcript|AAA47319.1|protein_coding|1/1|c.952A>C|p.Lys318Gln|952/3267|952/3267|318/1088||;BQB=0.572843;DP=36;DP4=19,7,3,5;GRANTHAM_SCORE=53;HOB=0.32;ICB=0.425;MQ=60;MQ0F=0;MQB=1;MQSB=1;RPB=0.658863;SGB=10.3229;VDB=0.693968	GT:PL	0/0:0,9,47	0/0:0,18,73	0/0:0,18,73	0/0:0,33,116	1/1:95,24,0
RF02	251	.	A	T	21.29	.	AC=2;AN=10;ANN=T|stop_gained|HIGH|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.235A>T|p.Lys79*|235/2643|235/2643|79/880||,T|upstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.-1371A>T|||||1371|WARNING_TRANSCRIPT_INCOMPLETE,T|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-1758A>T|||||1758|WARNING_TRANSCRIPT_NO_START_CODON;BQB=1;DP=24;DP4=18,0,6,0;GRANTHAM_SCORE=255;HOB=0.08;ICB=0.235294;MQ=60;MQ0F=0;MQB=1;RPB=0.566154;SGB=2.05141;VDB=0.0744703	GT:PL	0/0:0,15,57	0/1:31,0,5	0/1:31,0,5	0/0:0,9,42	0/0:0,24,69
```

END_DOC
 */
@Program(
	name="vcfgrantham",
	description="add grantham score from annotated VCF variant",
	keywords={"vcf","grantham"},
	creationDate = "20230503",
	modificationDate  = "20230503",
	menu="VCF Manipulation",
	jvarkit_amalgamion = true
	)
public class VcfGrantham extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.of(VcfGrantham.class);
	
	
	private abstract class Handler implements ToIntFunction<VariantContext>{
		abstract boolean isValid();
		}
	
	
	private class VepPredictionHandler extends Handler {
		private final VepPredictionParser parser;
		VepPredictionHandler(final VCFHeader header) {
			this.parser =  new VepPredictionParserFactory(header).get();
			}
		@Override boolean isValid() { return this.parser.isValid();}
		@Override
		public int applyAsInt(VariantContext ctx) {
			int score= 0;
			for(VepPredictionParser.VepPrediction pred : this.parser.getPredictions(ctx)) {
				final String s = pred.getAminoAcids();
				if(StringUtils.isBlank(s)) continue;
				int slash = s.indexOf("/");
				if(slash!=-1) {
					LOG.warning("cannot get '/' in  "+ pred.getAminoAcids());
					continue;
					}
				String aa1s = s.substring(0,slash);
				String aa2s = s.substring(slash+1);
				if(aa1s.length()!=1 || aa2s.length()!=1) {
					LOG.warning("cannot get 'x/y' in  "+ pred.getAminoAcids());
					continue;
					}
				if(aa1s.equals("*") || aa2s.equals("*")) {
					return GranthamScore.getDefaultScore();
					}
				
				int curr = GranthamScore.score(
						aa1s.charAt(0),
						aa2s.charAt(0)
						);
				score = Math.max(score, curr);
				
				}
			return score;
			}
		}
	
	private class BcsqPredictionHandler extends Handler {
		private final BcfToolsPredictionParser parser;
		BcsqPredictionHandler(final VCFHeader header) {
			this.parser =  new BcfToolsPredictionParserFactory(header).get();
			}
		@Override boolean isValid() { return this.parser.isValid();}
		@Override
		public int applyAsInt(VariantContext ctx) {
			int score = 0;
			for(BcfToolsPredictionParser.BcfToolsPrediction pred : this.parser.getPredictions(ctx)) {
				final String hgvs = pred.getAminoAcidChange();
				if(StringUtils.isBlank(hgvs)) continue;
				int lt = hgvs.indexOf(">");
				if(lt<0) {
					LOG.warning("cannot get '>' for "+hgvs);
					continue;
					}
				String aa1s = hgvs.substring(0,lt);
				String aa2s = hgvs.substring(lt+1);
				
				if(aa1s.endsWith("*") || aa2s.endsWith("*")) {
					return GranthamScore.getDefaultScore();
					}
				
				int curr = GranthamScore.score(
						aa1s.charAt(aa1s.length()-1),
						aa2s.charAt(aa2s.length()-1)
						);
				score = Math.max(score, curr);
				}
			return score;
			}
		}
	
	private class AnnPredictionHandler extends Handler {
		private final AnnPredictionParser parser;
		AnnPredictionHandler(final VCFHeader header) {
			this.parser =  new AnnPredictionParserFactory().header(header).get();
			}
		@Override boolean isValid() { return this.parser.isValid();}

		@Override
		public int applyAsInt(VariantContext ctx) {
			int score = 0;
			for(AnnPredictionParser.AnnPrediction pred : this.parser.getPredictions(ctx)) {
				final String hgvs = pred.getHGVSp();
				if(StringUtils.isBlank(hgvs)) continue;
				if(!hgvs.startsWith("p.")) continue;
				String s = hgvs.substring(2);
				// INDEL
				if(s.endsWith("fs")) return GranthamScore.getDefaultScore();
				
				
				if(s.startsWith("*") || s.endsWith("*")) {
					return GranthamScore.getDefaultScore();
					}
				
				if(s.endsWith("del")) {
					return GranthamScore.getDefaultScore();
					}
				
				if(s.endsWith("dup")) {
					return GranthamScore.getDefaultScore();
					}
				
				if(s.length()<3) {
					LOG.warning("small length for: "+hgvs);
					continue;
					}
				
				final String aa1s = s.substring(0,3);
				s=s.substring(3);
				final AminoAcid aa1 =AminoAcids.getAminoAcidFromThreeLettersCode(aa1s);
				if(aa1==null) {
					LOG.warning("cannot get aa1 for "+hgvs);
					continue;
					}
				final String aa2s = s.substring(s.length()-3);
				s=s.substring(0,s.length()-3);
				final AminoAcid aa2 =AminoAcids.getAminoAcidFromThreeLettersCode(aa2s);
				if(aa2==null) {
					LOG.warning("cannot get aa2 for "+hgvs);
					continue;
					}
				
				if(!StringUtils.isInteger(s)) {
					LOG.warning("cannot get loc for "+hgvs);
					continue;
					}
				final int curr = GranthamScore.score(aa1.getOneLetterCode(), aa2.getOneLetterCode());
				score = Math.max(score, curr);

				}
			return score;
			}
		}
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iter, VariantContextWriter out) {
		
		final VCFHeader header = iter.getHeader();
		final VCFInfoHeaderLine info = new VCFInfoHeaderLine(
				"GRANTHAM_SCORE",
				1,
				VCFHeaderLineType.Integer,
				"Max grantham score for this variant"
				);
		
		final List<Handler> handlers= Arrays.asList(
				new AnnPredictionHandler(header),
				new BcsqPredictionHandler(header),
				new VepPredictionHandler(header)
				).stream().
				filter(H->H.isValid()).
				collect(Collectors.toList());
	
		if(handlers.isEmpty()) {
			LOG.warn("no annotation handler was found for input file");
			}
		
		header.addMetaDataLine(info);
		JVarkitVersion.getInstance().addMetaData(this, header);
		out.writeHeader(header);
		while(iter.hasNext()) {
			final VariantContext ctx = iter.next();
			
			final int score = handlers.stream().
				mapToInt(H->H.applyAsInt(ctx)).
				max().
				orElse(0);
			
			if(score <=0 )
				{
				out.add(ctx);
				}
			else
				{
				out.add(new VariantContextBuilder(ctx).
						attribute(info.getID(), score).
						make());
				}
			
			}
		return 0;
		}
	
	public static void main(final String[] args) {
		new VcfGrantham().instanceMainWithExit(args);
	}

}
