package com.github.lindenb.jvarkit.tools.qqplot;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.OptionalDouble;
import java.util.function.ToDoubleFunction;

public class QQPlotter<T> {
	private List<PValueExtractor<T>> extractPValues;
	
	public interface PValueExtractor<U> {
		public String getLabel();
		public OptionalDouble apply(final U t);
	}
	
	private void doIt(Collection<T> data) {
		double min = data.stream().
				flatMap(T->extractPValues.stream().map(X->X.apply(T))).
				filter(T->T.isPresent()).
				mapToDouble(T->T.getAsDouble()).
				map(V-> -Math.log10(V)).
				min().orElse(1E-100);
		
		for(final PValueExtractor<T> extractor:this.extractPValues){
			double[] pobs = data.stream().
					map(T->extractor.apply(T)).
					filter(T->T.isPresent()).
					mapToDouble(T->T.getAsDouble()).
					map(V-> -Math.log10(V)).
					toArray();
			}
		}
	}
