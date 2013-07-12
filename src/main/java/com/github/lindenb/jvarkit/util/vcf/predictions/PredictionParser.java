package com.github.lindenb.jvarkit.util.vcf.predictions;

import java.util.List;

import org.broadinstitute.variant.variantcontext.VariantContext;

public interface PredictionParser {
public List<? extends Prediction> getPredictions(VariantContext ctx);
}
