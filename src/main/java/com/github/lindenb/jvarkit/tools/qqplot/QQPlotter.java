package com.github.lindenb.jvarkit.tools.qqplot;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.OptionalDouble;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

public class QQPlotter<T> {
	private List<PValueExtractor<T>> extractPValues;
	
	public interface PValueExtractor<U> {
		public String getLabel();
		public OptionalDouble apply(final U t);
	}
	private class DataPoint {
		public double getX() { return 0;}
		public double getY() { return 0;}
		}
	private class Series {
		private String label;
		private List<DataPoint> points;
		public String getLabel() { return "";}
		public List<DataPoint> getPoints() {return points;}
		}
	
	private void doIt(Collection<T> data) {
		final List<Series> series = new ArrayList<>(extractPValues.size());
		for(final PValueExtractor<T> extractor:this.extractPValues) {
			final Series serie = new Series();
			serie.label = extractor.getLabel();
			serie.points = data.stream().
					flatMap(T->extractPValues.stream().map(X->X.apply(T))).
					filter(T->T.isPresent()).
					mapToDouble(T->T.getAsDouble()).
					map(V-> -Math.log10(V)).
					mapToObj(V->new DataPoint()).
					collect(Collectors.toList());
			series.add(serie);
			}
		
		
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
