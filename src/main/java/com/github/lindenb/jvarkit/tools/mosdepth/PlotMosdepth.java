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
package com.github.lindenb.jvarkit.tools.mosdepth;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.RuntimeIOException;

/**
BEGIN_DOC

## Example
```
$ find . -type f > jeter.list
$ java -jar ~/src/jvarkit-git/dist/plotmosdepth.jar --max-coverage 100 --prefix 20210622.mosdepth. --format png jeter.list | R --vanilla > /dev/null
```

END_DOC
 */


@Program(name="plotmosdepth",
description="Plot Mosdepth output",
creationDate="20210621",
modificationDate="20210622",
keywords={"mosdepth"},
generate_doc=false
)
public class PlotMosdepth extends Launcher {
private static final Logger LOG = Logger.build(PlotMosdepth.class).make();
private final String SUFFIX_REGION_DIST = ".mosdepth.region.dist.txt";
private final String SUFFIX_GLOBAL_DIST = ".mosdepth.global.dist.txt";
private final String SUFFIX_SUMMARY = ".mosdepth.summary.txt";
private final String SUFFIX_REGION_BED_GZ = ".regions.bed.gz";

@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
private Path outputFile = null;
@Parameter(names={"--prefix"},description="output prefix.")
private String prefix="";
@Parameter(names={"--format"},description="output format.")
private Format outputFormat = Format.PDF;
@Parameter(names={"--max-coverage"},description="Max coverage fwhen plotting (percent of bases / coverage). Ignore if lower or equal to 0 ")
private int coverage_treshold = -1;
@Parameter(names={"--run-median"},description="runmed(coverage,'x') value for manhattan plot. Ignore if lower or equal to 0")
private int run_median_parameter = 5;
@Parameter(names={"--legend"},description="print legend")
private boolean print_legend=false;



private enum Where {GLOBAL,REGION};
private enum Format {PDF,PNG,SVG};

private static class RegionBed {
	String chrom;
	long pos;
	float cov;
	}

private static class RegionDist {
	String sample;
	String chrom;
	int cov;
	double frac;
	}
private static class Summary {
	String sample;
	double cov;
	
	}
private String quote(final String s) {
	return "\""+StringUtils.escapeC(s)+"\"";
	}

private String device(String fname,Map<String,String> hash) {
	String ext="";
	if(hash!=null && !hash.isEmpty()) {
		ext=","+hash.entrySet().stream().map(KV->KV.getKey()+"="+KV.getValue()).collect(Collectors.joining(","));
		}
	final String filename = this.prefix + fname;
	switch(outputFormat) {
		case SVG:
			return "svg("+quote(filename+".svg")+ext+")";
		case PNG:
			return "png("+quote(filename+".png")+ext+")";
		case PDF:
		default:
			return "pdf("+quote(filename+".pdf")+ext+")";
		}
	}


private void plotRegionsBed(final PrintWriter w,final Path path) {
		final String sample = getSample(path);
		try(BufferedReader br=IOUtils.openPathForBufferedReading(path)) {
				final List<RegionBed> points = br.lines().
				map(T->CharSplitter.TAB.split(T)).
				map(T->{
					final RegionBed bed = new RegionBed();
					bed.chrom = T[0];
					bed.pos  = Integer.parseInt(T[1]);
					bed.cov = Float.parseFloat(T[T.length-1]);
					return bed;
					}).
				collect(Collectors.toList());
			if(points.isEmpty()) return;
			final Set<String> chroms = new LinkedHashSet<>();
			RegionBed prev=null;
			long dx=0;
			for(int i=0;i< points.size();i++) {
				final RegionBed rb = points.get(i);
				chroms.add(rb.chrom);
				if(prev==null || !prev.chrom.equals(rb.chrom)) {
					dx = (prev==null?0:prev.pos);
					}
				rb.pos += dx;
				prev=rb;
				}
			for(final String c:chroms) {
				final double[] coverages = points.stream().
						filter(T->T.chrom.equals(c)).
						mapToDouble(X->X.cov).toArray();
				Arrays.sort(coverages);
				final int n= coverages.length;
				final double median;
				if(n%2==1)
					{ median = coverages[n/2]; } 
				else
					{ median = (coverages[(n-1)/2] + coverages[n/2])/2.0; }
				if(median==0) {
					points.removeIf(C->C.chrom.equals(c));
					continue;
					}
				points.stream().filter(T->T.chrom.equals(c)).forEach(T->T.cov/=(float)median);
				}
			
			final int max_cov = 3;
			Map<String,String> hash=new HashMap<String, String>();
			if(outputFormat.equals(Format.PNG)) hash.put("width", "2000");
			if(outputFormat.equals(Format.PDF)) hash.put("width", "21");
			w.println(device(sample+".manhattan",hash));
			w.println("max_cov_bed <- "+max_cov);
			int i=0;
			for(final String chrom: chroms) {
				w.print("X<-c(");
				w.print(points.stream().filter(P->P.chrom.equals(chrom)).map(T->String.valueOf(T.pos+".0")).collect(Collectors.joining(",")));
				w.println(")");
				
				w.print("Y<-c(");
				w.print(points.stream().filter(P->P.chrom.equals(chrom)).map(T->String.valueOf(T.cov)).collect(Collectors.joining(",")));
				w.println(")");
				if(this.run_median_parameter>0) {
					w.println("Y<-runmed(Y," + this.run_median_parameter + ")");
					}
				w.println("T2<-as.matrix(data.frame(X,Y))");
				w.println("head(T2)");

				String color = (i%2==0?"green":"red");
				if(chrom.equals("X") || chrom.equals("chrX")) {
					color="blue";
					}
				else if(chrom.equals("Y") || chrom.equals("chrY")) {
					color="pink";
					}
				
				final String xylim =
		                "xlim=c("+ points.stream().mapToLong(P->P.pos).min().orElse(0) + ","+
		                		points.stream().mapToLong(P->P.pos).max().orElse(1)+")," +
		                "ylim=c(0,max_cov_bed)";

				
				if(i==0) {
					w.println(
						"plot(T2,type = \"p\","
						+ "main=\""+sample+" Norm. Median depth\","+
		                "xlab=\"Genomic Position\","+
		                "ylab=\""+sample+" Norm Bed coverage\","+
		                "las=2,xaxt='n',"+
		                xylim+","+
		                "col= "+quote(color)+","+
		                "pch=3"+
		                ")");
					}
				else {
					w.println("par(new=TRUE)");
					w.println("plot(T2,type = \"p\",axes=FALSE,ann=FALSE,col=" +
						quote(color)+","+
						xylim+",pch=3)");
					}
				i++;
				}
			w.println("dev.off()");
			
			
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
		}
}

private void plotSummaries(PrintWriter w,final Where where,final List<Path> pathList) {
	if(pathList.isEmpty()) return;
	final List<Summary> points = new ArrayList<>(pathList.size());
	try {
		String percentile="average";
		for(final Path f:pathList) {
			final String sample = getSample(f);
			try(BufferedReader br=IOUtils.openPathForBufferedReading(f)) {
				String line;
				while((line=br.readLine())!=null) {
					final String tokens[] = CharSplitter.TAB.split(line);
					if(tokens[0].equals("chrom")) {
						percentile=tokens[3];
						continue;
						}
					if(tokens[0].equals(where.equals(Where.REGION)?"total_region":"total"))  {
						final Summary summary=new Summary();
						summary.sample=sample;
						summary.cov = Double.parseDouble(tokens[3]);
						//summary.min = Double.parseDouble(tokens[4]);
						//summary.max = Double.parseDouble(tokens[5]);
						points.add(summary);
						}
					}
				}
			}
		Collections.sort(points,(A,B)->Double.compare(A.cov, B.cov));
		if(points.size()<2) return;
		//w.println("par(old.par)");
		w.println(device("all.coverage"+(where.equals(Where.REGION)?".regions":""),Collections.emptyMap()));
		w.print("samples<-c(");
		w.print(points.stream().map(T->quote(T.sample)).collect(Collectors.joining(",")));
		w.println(")");
		w.print("cov<-c(");
		w.print(points.stream().map(T->String.valueOf(T.cov)).collect(Collectors.joining(",")));
		w.println(")");
		
		w.println(
				"barplot(cov,main="+quote(percentile + " Coverage"+(where.equals(Where.REGION)?" in selected regions":""))+","+
                "xlab="+quote("Sample")+","+
                "ylab="+quote("Coverage")+","+
                "las=2,"+
                "ylim=c(0,"+points.stream().mapToDouble(T->T.cov+1.0).max().orElse(1.0)+")," +
                "names.arg=samples"+
                ")");
		w.println("dev.off()");
		}
	catch(final Throwable err) {
		throw new RuntimeException(err);
		}
	}
private String getSample(final Path f) {
	final String sample = f.getFileName().toString();
	if(sample.endsWith(SUFFIX_REGION_DIST)) {
		return sample.substring(0,sample.length()-SUFFIX_REGION_DIST.length());
		}
	else if(sample.endsWith(SUFFIX_GLOBAL_DIST)) {
		return sample.substring(0,sample.length()-SUFFIX_GLOBAL_DIST.length());
		}
	else if( sample.endsWith(SUFFIX_SUMMARY)) {
		return sample.substring(0,sample.length()-SUFFIX_SUMMARY.length());
		}
	else if( sample.endsWith(SUFFIX_REGION_BED_GZ)) {
		return sample.substring(0,sample.length()-SUFFIX_REGION_BED_GZ.length());
		}
	else {
		return sample;
		}
	}
private List<RegionDist> readRegionDist(final Path f,final Predicate<RegionDist> accept) {
	final String sample = getSample(f);
	try(BufferedReader br=IOUtils.openPathForBufferedReading(f)) {
			return br.lines().
			map(T->CharSplitter.TAB.split(T)).
			map(T->{
				final RegionDist r = new RegionDist(); 
				r.chrom=T[0];
				r.sample = sample;
				r.cov = Integer.parseInt(T[1]);
				r.frac = Double.parseDouble(T[2]);
				if(!accept.test(r)) return null;
				return r;
				}).
			filter(T->T!=null).
			collect(Collectors.toList());
			}
		catch(final Throwable err) {
			throw new RuntimeException(err);
			}
		}

private void plotRegionDist(PrintWriter w,Where where ,final Path f) {
	final List<RegionDist> points = readRegionDist(f,X->!X.chrom.equals("total"));
	if(points.size()<1) return;
	final String sample =  getSample(f);
	final int max_cov = points.stream().filter(T->T.cov>0).mapToInt(T->T.cov).max().orElse(1);
	w.println(device(sample+(where.equals(Where.REGION)?".region":".global"),Collections.emptyMap()));
	final List<String> chromosomes = points.
			stream().
			map(T->T.chrom).
			collect(Collectors.toCollection(LinkedHashSet::new)).
			stream().
			collect(Collectors.toList());
	w.println("max_cov_dist <- " + (this.coverage_treshold>0?Math.min(100,max_cov):max_cov));
	w.println("colors <- rainbow("+chromosomes.size()+")");
	for(int i=0;i< chromosomes.size();i++) {
		final String chrom=chromosomes.get(i);
		w.print("cov<-c(");
		w.print(points.stream().filter(P->P.chrom.equals(chrom)).map(T->String.valueOf(T.cov)).collect(Collectors.joining(",")));
		w.println(")");
		w.print("frac<-c(");
		w.print(points.stream().filter(P->P.chrom.equals(chrom)).map(T->String.valueOf(T.frac)).collect(Collectors.joining(",")));
		w.println(")");
		w.println("T2<-as.matrix(data.frame(cov,frac))");
		if(i==0) {
			w.println(
				"plot(T2,type = \"l\","
				+ "main=\""+sample+(where.equals(Where.REGION)?" in selected region":"")+"\","+
                "xlab=\"Coverage\","+
                "ylab=\"Proportion of genome at coverage\","+
                "las=2,"+
                "xlim=c(0,max_cov_dist),"+
                "ylim=c(0,1.0),"+
                "col= colors[1],"+
                "pch=16"+
                ")");
			}
		else {
			w.println("par(new=TRUE)");
			w.print("plot(T2,type = \"l\",axes=FALSE,xlim=c(0,max_cov_dist),ylim=c(0,1.0),ann=FALSE,col=");
			w.print("colors["+(i+1)+"]");
                    	w.println(",pch=16)");
			}
		}
	if(print_legend) {
		w.println("legend(\"topright\",legend=c("+
			chromosomes.stream().map(S->quote(S)).collect(Collectors.joining(","))+
			"),title=\"chromosomes\",pch=16,col=colors)");
		}
	w.println("dev.off()");
	}

/** plot a BOX plot for each sample */
private void plotRegionsBed(PrintWriter w,final List<Path> paths) {
	if(paths.isEmpty()) return;
	final List<Map.Entry<String, double[]>> sample2covs = new ArrayList<>(paths.size());
	final Map<String,Double> sample2mean = new HashMap<>(paths.size());
	for(final Path path:paths) {
		try(BufferedReader br=IOUtils.openPathForBufferedReading(path)) {
			final double[] covs= br.lines().
				map(T->CharSplitter.TAB.split(T)).
				mapToDouble(T->Float.parseFloat(T[T.length-1])).
				toArray();
			if(covs.length==0) continue;
			final String sample = getSample(path);
			sample2covs.add(new AbstractMap.SimpleEntry<String,double[]>(sample,covs));
			sample2mean.put(sample, Arrays.stream(covs).average().orElse(0.0));
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	if(sample2covs.isEmpty()) return;
	Collections.sort(sample2covs,(A,B)->sample2mean.get(A.getKey()).compareTo(sample2mean.get(B.getKey())));
	
	w.println(device("boxplot.coverage",new HashMap<>()));
	w.println("boxplot(cov,");
	for(int i=0;i< sample2covs.size();i++) {
		if(i>0) w.print(",");
		double[] covs = sample2covs.get(i).getValue();
		w.print("c(");
		for(int j=0;j< covs.length;++j) {
			if(j>0) w.print(",");
			w.print(covs[j]);
			}
		w.println(")");
		}
	w.println(",main="+quote("Coverage in selected regions")+",ylab="+quote("Coverage"));
	w.print(",las=2,names=c("+ sample2covs.stream().map(KV->quote(KV.getKey())).collect(Collectors.joining(","))+ ")");
	w.println("dev.off()");
	}

private void plotTotalRegionDists(PrintWriter w,final List<Path> paths) {
	if(paths.isEmpty()) return;
	final boolean in_regions = paths.stream().anyMatch(F->F.getFileName().toString().endsWith(SUFFIX_REGION_DIST));
	final List<RegionDist> points =  
			paths.stream().
			flatMap(F->readRegionDist(F, (X->X.chrom.equals("total"))).stream()).
			collect(Collectors.toList())
			;
			
	final int max_cov = points.stream().filter(T->T.frac>0.0).mapToInt(T->T.cov).max().orElse(0) + 1;
	final List<String> samples = points.
			stream().
			map(T->T.sample).
			collect(Collectors.toCollection(TreeSet::new)).
			stream().
			collect(Collectors.toList());
	if(samples.size()<1) return;
	w.println(device(String.valueOf(samples.size())+"samples"+(in_regions?".region":".global"),Collections.emptyMap()));

	//w.println("par(old.par)");
	w.println("colors <- rainbow("+samples.size()+")");
	for(int i=0;i< samples.size();i++) {
		final String sample =samples.get(i);
		w.print("cov<-c(");
		w.print(points.stream().filter(P->P.sample.equals(sample)).map(T->String.valueOf(T.cov)).collect(Collectors.joining(",")));
		w.println(")");
		w.print("frac<-c(");
		w.print(points.stream().filter(P->P.sample.equals(sample)).map(T->String.valueOf(T.frac)).collect(Collectors.joining(",")));
		w.println(")");
		w.println("max_cov_dist <- " + (this.coverage_treshold>0?Math.min(100,max_cov):max_cov));
		
		final String xylim =   "xlim=c(0,max_cov_dist),ylim=c(0,1.0),";
		
		w.println("T2<-as.matrix(data.frame(cov,frac))");
		if(i==0) {
			w.println(
				"plot(T2,type = \"l\","
				+ "main=\"Coverage"+(in_regions?" in selected regions":"")+"\","+
                "xlab=\"Coverage\","+
                "ylab=\"Proportion of genome at coverage\","+
                "las=2,"+
                xylim +
                "col= colors[1],"+
                "pch=16"+
                ")");
			}
		else {
			w.println("par(new=TRUE)");
			w.print("plot(T2,type = \"l\","+xylim+"axes=FALSE,ann=FALSE,col=");
			w.print("colors["+(i+1)+"]");
                    	w.println(",pch=16)");
			}
		}
	
	if(print_legend) {
		w.println("legend(\"topright\",legend=c("+
				samples.stream().map(S->quote(S)).collect(Collectors.joining(","))+
			"),title=\"Sample\",pch=16,col=colors)");
		}
	w.println("dev.off()");
	}


@Override
public int doWork(final List<String> args) {
	try {
		final List<Path> files = IOUtils.unrollPaths(args);
		if(files.isEmpty()) {
			LOG.info("no file was provided");
			return -1;
			}
			
		
		try(PrintWriter pw = new PrintWriter(System.out)) {
			pw.println("# output to be piped in R");
			//see https://stackoverflow.com/questions/9292563
			//pw.println("old.par <- par(mar = c(0, 0, 0, 0))");
			
			for(int side=0;side<2;++side) {
				plotSummaries(pw,
						side==1?Where.GLOBAL:Where.REGION,
						files.stream().filter(F->F.getFileName().toString().
								endsWith(SUFFIX_SUMMARY)).
								collect(Collectors.toList())
						);
				}
		
			files.stream().
				filter(F->F.getFileName().toString().endsWith(SUFFIX_GLOBAL_DIST)).
				forEach(F->plotRegionDist(pw,Where.GLOBAL,F));

			
			files.stream().
					filter(F->F.getFileName().toString().endsWith(SUFFIX_REGION_DIST)).
					forEach(F->plotRegionDist(pw,Where.REGION,F));
			
			plotTotalRegionDists(pw,
					files.stream().
						filter(F->F.getFileName().toString().endsWith(SUFFIX_REGION_DIST)).
						collect(Collectors.toList())
					);
			plotTotalRegionDists(pw,
					files.stream().
						filter(F->F.getFileName().toString().endsWith(SUFFIX_GLOBAL_DIST)).
						collect(Collectors.toList())
					);
			files.stream().
				filter(F->F.getFileName().toString().endsWith(SUFFIX_REGION_BED_GZ)).
				forEach(F->plotRegionsBed(pw,F));

			
			plotRegionsBed(pw,
					files.stream().
						filter(F->F.getFileName().toString().endsWith(SUFFIX_GLOBAL_DIST)).
						collect(Collectors.toList())
					);
			
			pw.flush();
			}
		return 0;
		}
	catch(final Throwable err) {
		err.printStackTrace();
		return -1;
		}
	
	}


public static void main(final String[] args) {
		new PlotMosdepth().instanceMainWithExit(args);
		}
	}
