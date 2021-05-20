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
package com.github.lindenb.jvarkit.util.vcf.predictions;

import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
/**
    https://github.com/brentp/smoove/blob/master/annotate/annotate.go
 	vcf.AddInfoToHeader("smoove_gene", ".", "String", "genes overlapping variants. format is gene|feature:nfeatures:nbases,...")
 
 */
public class SmooveGenesParser {
	public static final String TAG="smoove_gene"; 
	private final boolean valid;

	public static interface Prediction {
		public String getGeneName();
		public String getFeature();
		public int getFeaturesCount();
		public int getBasesCount();
		public String getOriginalAttributeAsString();
		}

	private class PredictionImpl implements Prediction {
		String geneName;
		String feature;
		int featuresCount;
		int basesCount;
		@Override
		public String getGeneName() {
			return geneName;
			}
		@Override
		public String getFeature() {
			return feature;
			}
		@Override
		public int getFeaturesCount() {
			return featuresCount;
			}
		@Override
		public int getBasesCount() {
			return basesCount;
			}
		@Override
		public int hashCode() {
			int i= geneName.hashCode();
			i = i*31 + feature.hashCode();
			i= i*31 + Integer.hashCode(featuresCount);
			i= i*31 + Integer.hashCode(basesCount);
			return i;
			}
		@Override
		public boolean equals(Object obj) {
			if(obj==this)return true;
			if(obj==null || !(obj instanceof PredictionImpl)) return false;
			final PredictionImpl o = PredictionImpl.class.cast(obj);
			return this.geneName.equals(o.geneName) && 
					this.feature.equals(o.feature) && 
					this.featuresCount == o.featuresCount &&
					this.basesCount == o.basesCount
					;
			}
		@Override
		public String getOriginalAttributeAsString() {
			return this.geneName+"|"+this.feature+":"+featuresCount+":"+this.basesCount;
			}
		@Override
		public String toString() {
			return getOriginalAttributeAsString();
			}
		}


	public SmooveGenesParser(final VCFHeader header) {
		final VCFInfoHeaderLine h = header.getInfoHeaderLine(TAG);
		this.valid= h!=null;
		}

	public boolean isValid() {
		return valid;
		}

	public List<Prediction> parse(final VariantContext ctx) {
		if(!isValid() || ctx==null || !ctx.hasAttribute(getTag())) return Collections.emptyList();
		return ctx.getAttributeAsStringList(getTag(), null).
				stream().
				map(S->parseOne(S).orElse(null)).
				filter(S->S!=null).
				sorted((A,B)->A.getGeneName().compareTo(B.getGeneName())).
				collect(Collectors.toList())
				;
		}

	public Optional<Prediction> parseOne(final String s) {
		if(!isValid() || StringUtils.isBlank(s)) return Optional.empty();
		int pipe = s.indexOf('|');
		if(pipe==-1) return Optional.empty();
		int colon1 = s.indexOf(':',pipe+1);
		if(colon1==-1) return Optional.empty();
		int colon2 = s.indexOf(':',colon1+1);
		if(colon2==-1) return Optional.empty();

		final PredictionImpl p = new PredictionImpl();
		p.geneName = s.substring(0,pipe);
		p.feature = s.substring(pipe+1,colon1);
		p.featuresCount = Integer.parseInt(s.substring(colon1+1,colon2));
		p.basesCount = Integer.parseInt(s.substring(colon2+1));
		
		return Optional.of(p);
		}

	public String getTag() {
		return TAG;
		}
	}
