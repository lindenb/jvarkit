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
package com.github.lindenb.jvarkit.spliceai;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.CharSplitter;

import htsjdk.variant.variantcontext.VariantContext;

/** data extractor for spliceAI */
public interface SpliceAI {
static final String TAG= "SpliceAI";
public static String getTag() { return TAG;}

public String getAllele();
public String getGene();
public double getDeltaScoreAcceptorGain();
public double getDeltaScoreAcceptorLoss();
public double getDeltaScoreDonorGain();
public double getDeltaScoreDonorLoss();
public int getDeltaPositionAcceptorGain();
public int getDeltaPositionAcceptorLoss();
public int getDeltaPositionDonorGain();
public int getDeltaPositionDonorLoss();

public default double[] getDeltaScores() {
	return new double[] {
			getDeltaScoreAcceptorGain(),
			getDeltaScoreAcceptorLoss(),
			getDeltaScoreDonorGain(),
			getDeltaScoreDonorLoss()
		};
	}
public default double getDelatScoreMax() {
	return Arrays.stream(getDeltaScores()).max().orElse(0.0);
	}

public static List<SpliceAI> parse(final VariantContext vc) {
	if(!vc.hasAttribute(getTag())) return Collections.emptyList();
	return  vc.getAttributeAsStringList(getTag(), "").
			stream().
			map(S->CharSplitter.PIPE.split(S)).
			map(S->new SpliceAIImpl(S)).
			collect(Collectors.toList());
	}

/* ##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3.1 
 * variant annotation. These include delta scores (DS) and delta positions (DP) 
 * for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). 
 * Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">
 */
static class SpliceAIImpl implements SpliceAI {
	final String allele;
	final String gene;
	final double DS_AG,DS_AL,DS_DG,DS_DL;
	final int DP_AG,DP_AL,DP_DG,DP_DL;
	SpliceAIImpl(final String[] tokens) {
		this.allele = tokens[0];
		this.gene = tokens[1];
		this.DS_AG = Double.parseDouble(tokens[2]);
		this.DS_AL = Double.parseDouble(tokens[3]);
		this.DS_DG = Double.parseDouble(tokens[4]);
		this.DS_DL = Double.parseDouble(tokens[5]);

		this.DP_AG = Integer.parseInt(tokens[6]);
		this.DP_AL = Integer.parseInt(tokens[7]);
		this.DP_DG = Integer.parseInt(tokens[8]);
		this.DP_DL = Integer.parseInt(tokens[9]);
		}
	@Override
	public String getAllele() {
		return allele;
		}
	@Override public String getGene() {
		return gene;
		}
	@Override public double getDeltaScoreAcceptorGain() {
		return DS_AG;
		}
	@Override public double getDeltaScoreAcceptorLoss() {
		return DS_AL;
		}
	@Override public double getDeltaScoreDonorGain() {
		return DS_DG;
		}
	@Override public double getDeltaScoreDonorLoss() {
		return DS_DL;
		}
	
	@Override public int getDeltaPositionAcceptorGain() {
		return DP_AG;
		}
	@Override public int getDeltaPositionAcceptorLoss() {
		return DP_AL;
		}
	@Override public int getDeltaPositionDonorGain() {
		return DP_DG;
		}
	@Override public int getDeltaPositionDonorLoss() {
		return DP_DL;
		}
	
	@Override
	public String toString() {
		return allele+"|"+gene+"|"+DS_AG+"|"+DS_AL+"|"+DS_DG+"|"+DS_DL+
				"|"+DP_AG+"|"+DP_AL+"|"+DP_DG+"|"+DP_DL;
		}
	}
}
