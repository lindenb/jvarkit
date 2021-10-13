package com.github.lindenb.jvarkit.tools.bcftools;

import java.io.BufferedReader;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;
import java.util.function.ToLongFunction;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.RangeOfDoubles;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
BEGIN_DOC

## Example
```
$ find . -type f > jeter.list
$ java -jar ~/src/jvarkit-git/dist/plotmosdepth.jar --max-coverage 100 --prefix 20210622.mosdepth. --format png jeter.list | R --vanilla > /dev/null
```

END_DOC
 */


@Program(name="plotbcftoolsstats",
description="Plot bcftools stats output",
creationDate="20210622",
modificationDate="20210622",
keywords={"bcftools","qc","vcf"},
generate_doc=false
)
public class PlotBcftoolsStats extends Launcher {
	private static final Logger LOG = Logger.build(PlotBcftoolsStats.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--prefix"},description="output prefix.")
	private String prefix="";
	@Parameter(names={"--format"},description="output format.")
	private Format outputFormat = Format.PDF;
	@Parameter(names={"--categories","--phenotypes"},description="A tab delimited file (sample)<tab>(category)")
	private Path sample2catPath = null;
	@Parameter(names={"--sections"},description="Limit to those sections. eg. 'QUAL,PSC,PSI'. Comma-separated values. Ignore if empty. Inverse the selection if it starts with '^'.")
	private String sectionsSelectStr = "";

	@DynamicParameter(names = "-D", description = "set some css style elements. '-Dkey=value'. Undocumented.")
	private Map<String, String> dynaParams = new HashMap<>();


	private enum Format {PDF,PNG,SVG};
	private final Map<String,String> id2vcfs = new TreeMap<>();
	private final Map<String, List<String>> phenotype2samples = new TreeMap<>();
	private final Map<String,String> sample2phenotype = new TreeMap<>();
	private int depth_bin_size = -1;

	private abstract class Section {
		final String ID;
		final String header;
		final int expectNCols;
		boolean enabled =  true;
		protected Section(String ID,String header) {
			this.ID = ID;
			this.header = header;
			this.expectNCols= CharSplitter.TAB.split(this.header).length;
			}
		final List<List<String>> lines = new ArrayList<>();
		
		Set<String> getIds() {
			return this.lines.stream().map(L->L.get(1)).collect(Collectors.toCollection(LinkedHashSet::new));
			}
		
		List<List<String>> getLinesForId(final String id) {
			return this.lines.stream().
					filter(L->L.get(1).equals(id)).
					collect(Collectors.toList());
			}
		
		void add(final String[] tokens) {
			if(tokens.length!=this.expectNCols) {
				throw new IllegalArgumentException("boum");
				}
			lines.add(Arrays.asList(tokens));
			}
		abstract void print(PrintWriter w);
		}
	
	/** SiS, Singleton stats *******************************************************************/
	private class SiSSection extends Section {
		SiSSection() {
			super("SiS","# SiS\t[2]id\t[3]allele count\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable");
			}
		@Override
		void print(PrintWriter w) {
		}
		}
	
	/** TSTV transitions/transversions *******************************************************************/
	private class TSTVSection extends Section {
		TSTVSection() {
			super("TSTV","# TSTV\t[2]id\t[3]ts\t[4]tv\t[5]ts/tv\t[6]ts (1st ALT)\t[7]tv (1st ALT)\t[8]ts/tv (1st ALT)");
			}
		@Override
		void print(PrintWriter w) {
			for(String id : getIds()) {
				final List<List<String>> lines =getLinesForId(id);
				w.println(device(this.ID+"."+id,null));
				final double ts = Double.parseDouble(lines.get(0).get(2));
				final double tv = Double.parseDouble(lines.get(0).get(3));
				final double ts1 = Double.parseDouble(lines.get(0).get(5));
				final double tv1 = Double.parseDouble(lines.get(0).get(6));

				w.println(
						"barplot(c("+ts+","+tv+","+ts1+","+tv1+")," +
						"main="+quote("TS/TV")+","+
						"sub="+quote(id2vcfs.getOrDefault(id,""))+","+
		                "ylab=\"Count\","+
		                "xlab=\"\","+
		                "ylim=c(0,"+Math.max(ts,tv)+"),"+
		                "names=c("+quote("TS")+","+quote("TV")+","+quote("TS (1st ALT)")+","+quote("TV (1st ALT)")+"),"+
		                "col=c("+quote("green")+","+quote("blue")+","+quote("green")+","+quote("blue")+"),"+
						"las=0)");
				
				w.println("dev.off()");
				}
			
			}
		}
	/** ICS Indel context summary *******************************************************************/
	private class ICSSection extends Section {
		ICSSection() {
			super("ICS","# ICS\t[2]id\t[3]repeat-consistent\t[4]repeat-inconsistent\t[5]not applicable\t[6]c/(c+i) ratio");
			}
		@Override
		void print(PrintWriter w) {
			for(String id : getIds()) {
				final List<List<String>> lines =getLinesForId(id);
				w.println(device(this.ID+"."+id,null));
				final long repeat_consistent = Long.parseLong(lines.get(0).get(2));
				final long repeat_inconsistent = Long.parseLong(lines.get(0).get(3));
				final long no_applicable = Long.parseLong(lines.get(0).get(4));

				w.println(
						"barplot(c("+repeat_consistent+","+repeat_inconsistent+","+no_applicable+")," +
						"main="+quote("Indel context summary")+","+
						"sub="+quote(id2vcfs.getOrDefault(id,""))+","+
		                "ylab=\"Count\","+
		                "xlab=\"\","+
		                "ylim=c(0,"+Math.max(no_applicable,Math.max(repeat_consistent,repeat_inconsistent))+"),"+
		                "names=c("+quote("Consistent")+","+quote("Inconsistent")+","+quote("Not Applicable")+"),"+
						"las=0)");
				
				w.println("dev.off()");
				}
			
			}
		}

	
	/** ICL  Indel context by length *******************************************************************/
	private class ICLSection extends Section {
		ICLSection() {
			super("ICL","# ICL\t[2]id\t[3]length of repeat element\t[4]repeat-consistent deletions)\t[5]repeat-inconsistent deletions\t[6]consistent insertions\t[7]inconsistent insertions\t[8]c/(c+i) ratio");
			}
		@Override
		void print(PrintWriter w) {
			for(String id : getIds()) {
				final List<List<String>> lines =getLinesForId(id);
				w.println(device(this.ID+"."+id,null));
				for(int i=0;i< lines.size();i++) {
					w.println("r"+i+" <- c("+ lines.get(i).subList(3, 8).stream().collect(Collectors.joining(",")) +")");
					}
				
				w.println("T2 <- as.matrix(data.frame("+ IntStream.range(0, lines.size()).mapToObj(i->"r"+i).collect(Collectors.joining(",")) +"))");
				w.println(
						"barplot(T2," +
						"main="+quote("Indel context summary")+","+
						"sub="+quote(id2vcfs.getOrDefault(id,""))+","+
		                "ylab=\"Count\","+
		                "xlab=\"length of repeat element\","+
		                "names=c("+lines.stream().map(L->quote(L.get(2))).collect(Collectors.joining(","))+"),"+
		                "legend=c("+String.join(",", quote("repeat-consistent deletions"),quote("repeat-inconsistent deletions"),quote("repeat-consistent insertions"),quote("repeat-inconsistent insertions"))+"),"+
		                "beside=TRUE,"+
		                "col=rainbow(4),"+
						"las=0)");
				
				w.println("dev.off()");
				}
			}
		}
	
	/** IDD InDel distribution *******************************************************************/
	private class IDDSection extends Section {
		IDDSection() {
			super("IDD","# IDD\t[2]id\t[3]length (deletions negative)\t[4]number of sites\t[5]number of genotypes\t[6]mean VAF");
			}
		@Override
		void print(PrintWriter w) {
			print(w,"Number of Sites",L->Long.parseLong(L.get(3)));
			print(w,"Number of Genotypes",L->Long.parseLong(L.get(4)));
			print(w,"Mean VAF",L->Double.parseDouble(L.get(5)));
			}
		
		
		private void  print(PrintWriter w,String title,Function<List<String>,Number> extractor) {
			final String fname= title.replaceAll(" of ", "").replaceAll("\\s+","");
			for(String id : getIds()) {
				final List<List<String>> lines =getLinesForId(id);
				w.println(device(this.ID+"."+fname +"."+id,null));
				w.println(
						"plot("+
						"x=c("+ lines.stream().map(L->L.get(2)).collect(Collectors.joining(","))+"),"+
						"y=c("+ lines.stream().map(L->String.valueOf(extractor.apply(L))).collect(Collectors.joining(","))+")," +
						"main="+quote("InDel distribution: "+ title)+","+
						"sub="+quote(id2vcfs.getOrDefault(id,""))+","+
						"type=\"l\","+
		                "xlab=\"length\","+
		                "ylab="+quote("Count "+title)+
		                ")"
						);
						
				w.println("dev.off()");
				}
			}
		}

	
	/** AF Stats by non-reference allele frequency *******************************************************************/
	private class AFSection extends Section {
		AFSection() {
			super("AF","# AF\t[2]id\t[3]allele frequency\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable");
			}
		@Override
		void print(PrintWriter w) {
			print(w,"Number of SNPs",L->Long.parseLong(L.get(3)));
			print(w,"Number of Transitions",L->Long.parseLong(L.get(4)));
			print(w,"Number of Transversions",L->Long.parseLong(L.get(5)));
			print(w,"Number of Indels",L->Long.parseLong(L.get(6)));
			print(w,"Repeat Consistent",L->Long.parseLong(L.get(7)));
			print(w,"Repeat Inconsistent",L->Long.parseLong(L.get(8)));
			print(w,"Not applicable",L->Long.parseLong(L.get(8)));
			}
		
		
		private void  print(PrintWriter w,String title,ToLongFunction<List<String>> extractor) {
			final String fname= title.replaceAll(" of ", "").replaceAll("\\s+","");
			for(String id : getIds()) {
				final List<List<String>> lines =getLinesForId(id);
				w.println(device(this.ID+"."+fname +"."+id,null));
				w.println(
						"barplot(c("+
						lines.stream().map(L->String.valueOf(extractor.applyAsLong(L))).collect(Collectors.joining(",")) +
						"),main="+quote("Stats by non-reference allele frequency: "+ title)+","+
						"sub="+quote(id2vcfs.getOrDefault(id,""))+","+
		                "xlab=\"AF\","+
		                "ylab="+quote("Count "+title)+","+
		                "ylim=c(0,"+lines.stream().mapToLong(extractor).max().orElse(1L)+"),"+
		                "names=c("+lines.stream().map(L->L.get(2)).collect(Collectors.joining(","))+"),"+
						"las=2)");
						
				w.println("dev.off()");
				}
			}
		}
	/** HWE Section *******************************************************************/
	private class HWESection extends Section {
		HWESection() {
			super("HWE","# HWE\t[2]id\t[3]1st ALT allele frequency\t[4]Number of observations\t[5]25th percentile\t[6]median\t[7]75th percentile");
			}
		@Override
		void print(PrintWriter w) {
			print(w,"Number of observations",L->Long.parseLong(L.get(3)));
			print(w,"25th Percentile",L->Double.parseDouble(L.get(4)));
			print(w,"Median",L->Double.parseDouble(L.get(5)));
			print(w,"75th Percentile",L->Double.parseDouble(L.get(6)));
			}
		
		
		private void  print(PrintWriter w,String title,Function<List<String>,Number> extractor) {
			final String fname= title.replaceAll(" of ", "").replaceAll("\\s+","");
			for(String id : getIds()) {
				final List<List<String>> lines =getLinesForId(id);
				w.println(device(this.ID+"."+fname +"."+id,null));
				w.println(
						"barplot(c("+
						lines.stream().map(L->String.valueOf(extractor.apply(L))).collect(Collectors.joining(",")) +
						"),main="+quote("HWE: "+ title)+","+
						"sub="+quote(id2vcfs.getOrDefault(id,""))+","+
		                "xlab=\"AF\","+
		                "ylab="+quote(title)+","+
		                "ylim=c(0,"+lines.stream().map(extractor).mapToDouble(X->X.doubleValue()).max().orElse(1L)+"),"+
		                "names=c("+lines.stream().map(L->L.get(2)).collect(Collectors.joining(","))+"),"+
						"las=2)");
						
				w.println("dev.off()");
				}
			}
		}
	/** Depth distribution *******************************************************************/
	private class DPSection extends Section {
		DPSection() {
			super("DP","# DP\t[2]id\t[3]bin\t[4]number of genotypes\t[5]fraction of genotypes (%)\t[6]number of sites\t[7]fraction of sites (%)");
			}
		@Override
		void print(PrintWriter w) {
			if(depth_bin_size<1) {
				LOG.warn("I was not able to find the --depth bin-size");
				return;
				}
			print(w,"Number of Genotypes",L->Double.parseDouble(L.get(3)));
			print(w,"Fraction of Genotypes",L->Double.parseDouble(L.get(4)));
			print(w,"Number of Sites",L->Double.parseDouble(L.get(5)));
			print(w,"Fraction of Sites",L->Double.parseDouble(L.get(6)));
			}
		
		private String binToLabel(String s) {
			if(s.startsWith(">")) return quote(s);
			int n = Integer.parseInt(s);
			return quote("["+(n*depth_bin_size)+"-"+((n+1)*depth_bin_size)+"[");
		}
		
		private void  print(PrintWriter w,String title,ToDoubleFunction<List<String>> extractor) {
			final String fname= title.replaceAll(" of ", "").replaceAll("\\s+","");
			for(String id : getIds()) {
				final List<List<String>> lines =getLinesForId(id);
				w.println(device(this.ID+"."+fname +"."+id,null));
				w.println(
						"barplot(c("+
						lines.stream().map(L->String.valueOf(extractor.applyAsDouble(L))).collect(Collectors.joining(",")) +
						"),main="+quote("DEPTH: "+ title)+","+
						"sub="+quote(id2vcfs.getOrDefault(id,""))+","+
		                "xlab=\"Depth\","+
		                "ylab="+quote(title)+","+
		                "ylim=c(0,"+lines.stream().mapToDouble(extractor).max().orElse(1.0)+"),"+
		                "names=c("+lines.stream().map(L->binToLabel(L.get(2))).collect(Collectors.joining(","))+"),"+
						"las=2)");
						
				w.println("dev.off()");
				}
			}
		}
	/** ST (Substitution types) *******************************************************************/
	private class STSection extends Section {
		STSection() {
			super("ST","# ST\t[2]id\t[3]type\t[4]count");
			}
		
		private String substitutionToColor(final String s) {
			if(s.length()!=3) throw new IllegalArgumentException("not [ATGC]>[ATGC]");
			String color ="red";
			char b1 = s.charAt(0);
			char b2 = s.charAt(2);
			if(AcidNucleics.isPuryn(b1) && AcidNucleics.isPuryn(b2)) {
				color = "blue";
			} else if(AcidNucleics.isPyrimidic(b1) && AcidNucleics.isPyrimidic(b2)) {
				color = "green";
			} 
			
			return quote(color);
		}
		
		@Override
		void print(PrintWriter w) {
			for(String id : getIds()) {
				w.println(device(this.ID+"."+id,null));
				final Counter<String> counter = new Counter<String>();
				for(final List<String> line:getLinesForId(id)) {
					if(line.get(3).equals("0")) continue;
					counter.incr(line.get(2),Long.parseLong(line.get(3)));
					}
				if(counter.isEmpty()) continue;

				w.println(
						"barplot(c("+
						counter.keySetDecreasing().stream().map(K->String.valueOf(counter.count(K))).collect(Collectors.joining(",")) +
						"),main=\"Substitution types\","+
						"sub="+quote(id2vcfs.getOrDefault(id,""))+","+
		                "xlab=\"Subsitution\","+
		                "ylab=\"Count Variants\","+
		                "ylim=c(0,"+counter.getMaxCount().orElse(1L)+"),"+
		                "names=c("+counter.keySetDecreasing().stream().map(S->quote(S)).collect(Collectors.joining(","))+"),"+
						"col= c(" + 
						counter.keySetDecreasing().stream().map(K->substitutionToColor(K)).collect(Collectors.joining(",")) +
						"),las=2)");
						}
			w.println("legend(\"topright\","+ 
				       "legend = c(\"purin-purin\",\"pyr-pyr\",\"other\"),"+ 
				       "fill = c("+substitutionToColor("A>G")+","+substitutionToColor("T>C")+","+substitutionToColor("A>C")+"))");
				
			w.println("dev.off()");
			}
		}
	
	
	/** QUAL *******************************************************************/
	private class QualSection extends Section {
		QualSection() {
			super("QUAL","# QUAL\t[2]id\t[3]Quality\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)\t[7]number of indels");
			}
		@Override
		void print(PrintWriter w) {
			final RangeOfDoubles ranges = RangeOfDoubles.fromTo(0,
					Double.parseDouble(dynaParams.getOrDefault("qual.max","1000")),
					Double.parseDouble(dynaParams.getOrDefault("qual.step","100"))
					); 
			for(String id : getIds()) {
				w.println(device(this.ID+"."+id,null));
				List<List<String>> L = getLinesForId(id);
				long max_count=1L;
				List<String> colNames= new ArrayList<>();
				for(RangeOfDoubles.Range r: ranges.getRanges()) {
					final String colName = "count" + colNames.size();
					boolean first= true;
					w.print(colName+"<-c(");
					for(int i=0;i< 4;i++) {
						long n = 0L;
						for(final List<String> line:L) {
							final double qual = Double.parseDouble(line.get(2));
							if(!r.contains(qual)) continue;
							n += Long.parseLong(line.get(3+i));
							}
						if(!first) w.print(",");
						first=false;
						w.print(n);
						max_count = Math.max(n,max_count);
						}
					w.println(")");
					colNames.add(colName);
					}
				w.println("T2 <- as.matrix(data.frame("+String.join(",", colNames) +"))");
				
				device(this.ID+"."+id,null);
				w.println(
					"barplot(T2,main=\"Stats by quality\","+
					"sub="+quote(id2vcfs.getOrDefault(id,""))+","+
	                "xlab=\"Quality\","+
	                "ylab=\"Count Variants\","+
	                "ylim=c(0,"+max_count+"),"+
	                "names=c("+ranges.stream().map(S->quote(S.toString())).collect(Collectors.joining(","))+"),"+
	                "legend=c("+Arrays.stream(CharSplitter.TAB.split(this.header)).skip(3).map(S->quote(removeColumn(S))).collect(Collectors.joining(","))+
					"),col= rainbow(4),beside=TRUE,las=2)");
					}

				w.println("dev.off()");
				}
			}
	/** SN *******************************************************************/
	private class SNSection extends Section {
		SNSection() {
			super("SN","# SN\t[2]id\t[3]key\t[4]value");
			}
		@Override
		void print(PrintWriter w) {
			for(String id : getIds()) {
				w.println(device(this.ID+"."+id,null));
				final List<List<String>> L = getLinesForId(id).stream().
						filter(S->!S.get(2).contains("number of samples")).
						collect(Collectors.toList());
	
				
				device(this.ID+"."+id,null);
				w.println(
					"barplot(c("+L.stream().map(T->String.valueOf(T.get(3))).collect(Collectors.joining(","))+")," +
					"main=\"Summary Numbers\","+
					"sub="+quote(id2vcfs.getOrDefault(id,""))+","+
	                "ylab=\"Count Variants\","+
	                "names=c("+L.stream().map(T->quote(removeColumn(T.get(2)))).collect(Collectors.joining(","))+"),"+
					"las=2,cex.names=0.5)");
					}

				w.println("dev.off()");
				}
			}
	
	private abstract class PerSampleSection extends Section {
		protected class Column {
			final String label;
			final ToDoubleFunction<List<String>>  fun;
			Column(final String label,ToDoubleFunction<List<String>>  fun) {
				this.label= label;
				this.fun = fun;
				}
			}
		PerSampleSection(String id,String header) {
			super(id,header);
			}
		
		protected void groupByPhenotype(PrintWriter w,String id,String title,Function<List<String>, Number> fun) {
			if(sample2phenotype.isEmpty()) return;
			final List<List<String>> L = getLinesForId(id);

			final Set<String> samples = L.stream().map(S->S.get(2)).collect(Collectors.toCollection(HashSet::new));
			samples.retainAll(sample2phenotype.keySet());
			
			final List<String> phenotypes = new ArrayList<>(samples.stream().map(S->sample2phenotype.get(S)).
					collect(Collectors.toCollection(TreeSet::new)));
			if(phenotypes.isEmpty()) {
				LOG.info("no phenotype for "+samples);
				return;
			}
			
			w.println(device(this.ID+"."+id+".boxplot."+title.toLowerCase().replace(' ', '_'),null));
			for(int i=0;i< phenotypes.size();i++) {
				final String phenotype = phenotypes.get(i);
				w.println("p"+i+" <- c("+
					L.stream().
						filter(X->phenotype.equals(sample2phenotype.get(X.get(2)))).
						map(X->String.valueOf(fun.apply(X))).
						collect(Collectors.joining(",")) +
					")");
				}
			w.println("boxplot("+
				IntStream.range(0, phenotypes.size()).mapToObj(X->"p"+X).collect(Collectors.joining(","))+","+
				"main="+quote(title)+",ylab="+quote(title)+","+
				"las=2,"+
				"names=c("+ phenotypes.stream().map(S->quote(S)).collect(Collectors.joining(","))+")"+
				")");
			w.println("dev.off()");
			}
		
		protected void print(PrintWriter w,String id,String title,String extra,final List<Column> columns) {
			final boolean has_phenotype = columns.size()==1 && !sample2phenotype.isEmpty();

			final List<List<String>> L;
			
			if(columns.size()==1) {
				L = getLinesForId(id).stream().sorted((A,B)->Double.compare(columns.get(0).fun.applyAsDouble(B), columns.get(0).fun.applyAsDouble(A))).collect(Collectors.toList());
				}
			else
				{
				L = getLinesForId(id);
				}
			
			if(columns.stream().flatMapToDouble(C->L.stream().mapToDouble(C.fun)).allMatch(D->D==0)) return;
			
			if(has_phenotype) {
				w.println("pc <- rainbow("+phenotype2samples.size()+")");
				}
			
			final UnaryOperator<String> pheno2color= PH ->{
				if(!StringUtils.isBlank(PH)) {
					int i=0;
					for(String sp:phenotype2samples.keySet()) {
						if(sp.equals(PH)) break;
						i++;
						}
					if(i==phenotype2samples.size()) throw new IllegalArgumentException("Cannot find "+PH+" in "+phenotype2samples.keySet());
					return "pc["+(i+1)+"]";
					}
				return  "\"gray\"";
				};
			
			final UnaryOperator<String> sample2color= S -> pheno2color.apply(sample2phenotype.get(S));
			
			final List<String> colNames = new ArrayList<>();
			final HashMap<String, String> hash = new HashMap<>();
			hash.put("width", String.valueOf((int)Math.max(7.0,L.size()*7/40.0)));
			w.println(device(this.ID+"."+id+"."+title.toLowerCase().replace(' ', '_'),hash));
			//loop over samples
			for(int i=0;i< L.size();i++) {
				final List<String> row= L.get(i);
				final String colName = "c"+colNames.size();
				w.println("# sample "+row.get(2));
				w.print(colName+"<-c(");
				w.print(columns.stream().map(C->String.valueOf(C.fun.applyAsDouble(row))).collect(Collectors.joining(",")));
				w.println(")");
				colNames.add(colName);
				}
			
			if(has_phenotype) {
				w.println("T2 <- c("+String.join(",", colNames) +")");
				}
			else {
				w.println("T2 <- as.matrix(data.frame("+String.join(",", colNames) +"))");
				}

			w.println(
					"barplot(T2,main="+quote(title)+","+
					"sub="+quote(id2vcfs.getOrDefault(id,""))+","+
	                "xlab=\"Sample\","+
	                "ylab=\"Count Variants\","+
	                "names=c("+L.stream().map(S->quote(S.get(2))).collect(Collectors.joining(","))+"),"+
	                (has_phenotype?"":"legend=c("+columns.stream().map(S->quote(S.label)).collect(Collectors.joining(","))+"),")+
	                "col="+(has_phenotype?"c("+ L.stream().map(S->sample2color.apply(S.get(2))).collect(Collectors.joining(",")) +")":"rainbow("+ columns.size() +")")+","+
					"las=2"+extra+")");
			
			if(has_phenotype) {
				w.println("legend(\"topright\","+ 
					"legend =c("+ phenotype2samples.keySet().stream().map(S->quote(S)).collect(Collectors.joining(","))+"),"+ 
					"fill = c("+phenotype2samples.keySet().stream().map(S->pheno2color.apply(S)).collect(Collectors.joining(","))+")"+
					")");
				}
			
			w.println("dev.off()");
			}

	}
	
	/** PSC per sample count *******************************************************************/
	private class PSCSection extends PerSampleSection {
		PSCSection() {
			super("PSC","# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions\t[8]nTransversions\t[9]nIndels\t[10]average depth\t[11]nSingletons\t[12]nHapRef\t[13]nHapAlt\t[14]nMissing");
			}
		@Override
		void print(PrintWriter w) {
			for(String id : getIds()) {
				print(w,id,"Average Depth","",Arrays.asList(new Column("Depth",L->Double.parseDouble(L.get(9)))));
				print(w,id,"Missing","",Arrays.asList(new Column("Missing",L->Double.parseDouble(L.get(13)))));
				print(w,id,"Singletons","",Arrays.asList(new Column("Singletons",L->Double.parseDouble(L.get(10)))));
				print(w,id,"Genotypes","",
						Arrays.asList(
								new Column("HOM_REF",L->Double.parseDouble(L.get(3))),
								new Column("HET",L->Double.parseDouble(L.get(5))),
								new Column("HOM_VAR",L->Double.parseDouble(L.get(4))),
								new Column("NO_CALL",L->Double.parseDouble(L.get(13)))
								));
				
				
				groupByPhenotype(w,id,"Average Depth",L->Double.parseDouble(L.get(9)));
				groupByPhenotype(w,id,"Missing",L->Double.parseDouble(L.get(13)));
				groupByPhenotype(w,id,"Het",L->Double.parseDouble(L.get(5)));
				groupByPhenotype(w,id,"HOM_VAR",L->Double.parseDouble(L.get(4)));
				groupByPhenotype(w,id,"Singletons",L->Double.parseDouble(L.get(10)));
				}
			}
		}
	
	/** PSC per sample Indels *******************************************************************/
	private class PSISection extends PerSampleSection {		
		
		PSISection() {
			super("PSI","# PSI\t[2]id\t[3]sample\t[4]in-frame\t[5]out-frame\t[6]not applicable\t[7]out/(in+out) ratio\t[8]nInsHets\t[9]nDelHets\t[10]nInsAltHoms\t[11]nDelAltHoms");
			}
		@Override
		void print(PrintWriter w) {
			for(String id : getIds()) {
				print(w,id,"In Frame","",Arrays.asList(new Column("In Frame",L->Double.parseDouble(L.get(3)))));
				print(w,id,"Out Frame","",Arrays.asList(new Column("Out Frame",L->Double.parseDouble(L.get(4)))));
				print(w,id,"Not applicable","",Arrays.asList(new Column("Not applicable",L->Double.parseDouble(L.get(5)))));
				print(w,id,"Out Fraction","",Arrays.asList(new Column("Out Fraction",L->Double.parseDouble(L.get(6)))));
				print(w,id,"nInsHets","",Arrays.asList(new Column("nInsHets",L->Double.parseDouble(L.get(7)))));
				print(w,id,"nDelsHets","",Arrays.asList(new Column("nDelsHets",L->Double.parseDouble(L.get(8)))));
				print(w,id,"nInsAltHoms","",Arrays.asList(new Column("nInsAltHoms",L->Double.parseDouble(L.get(9)))));
				print(w,id,"nDelAltHoms","",Arrays.asList(new Column("nDelAltHoms",L->Double.parseDouble(L.get(10)))));
				
				groupByPhenotype(w,id,"nInsHets",L->Double.parseDouble(L.get(7)));
				groupByPhenotype(w,id,"nDelsHets",L->Double.parseDouble(L.get(8)));
				groupByPhenotype(w,id,"nInsAltHoms",L->Double.parseDouble(L.get(9)));
				groupByPhenotype(w,id,"nDelAltHoms",L->Double.parseDouble(L.get(10)));
				}
			}
		}
	
	private String quote(final String s) {
		return "\""+StringUtils.escapeC(s)+"\"";
		}

	private String removeColumn(final String s) {
		if(!s.startsWith("[")) return s;
		int i=s.indexOf("]");
		if(i<=0) return s;
		return s.substring(i+1);
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


	@Override
	public int doWork(final List<String> args) {
		try {
			final List<Section> sections = Arrays.asList(
				new AFSection(),
				new DPSection(),
				new HWESection(),
				new ICSSection(),
				new ICLSection(),
				new IDDSection(),
				new STSection(),
				new QualSection(),
				new SNSection(),
				new PSCSection(),
				new PSISection(),
				new TSTVSection(),
				new SiSSection()
				);
			
			if(!StringUtils.isBlank(this.sectionsSelectStr)) {
				final boolean inverse = this.sectionsSelectStr.startsWith("^");
				final Set<String> ids = Arrays.stream(CharSplitter.COMMA.split(inverse?this.sectionsSelectStr.substring(1):this.sectionsSelectStr)).
					map(S->S.trim()).
					filter(S->!StringUtils.isBlank(S)).
					collect(Collectors.toSet());
				for(Section sec: sections) {
					if(inverse) {
						sec.enabled = !ids.contains(sec.ID);
						}
					else
						{
						sec.enabled = ids.contains(sec.ID);
						}
					}
				}
			
			if(sample2catPath!=null) {
				try(BufferedReader br= IOUtils.openPathForBufferedReading(sample2catPath)) {
					String line;
					while((line=br.readLine())!=null) {
						
						if(line.startsWith("#") || StringUtils.isBlank(line)) {
							continue;
							}
						final String[] tokens = CharSplitter.TAB.split(line);
						if(tokens.length!=2) throw new JvarkitException.TokenErrors(2, tokens);
						if(sample2phenotype.containsKey(tokens[0])) {
							LOG.error("multiple phenotypes for "+tokens[0]+" in "+sample2catPath);
							return -1;
							}
						sample2phenotype.put(tokens[0], tokens[1]);
						List<String> samples = phenotype2samples.get(tokens[1]);
						if(samples==null) {
							samples = new ArrayList<>();
							phenotype2samples.put(tokens[1],samples);
							}
						samples.add(tokens[0]);
					}
				}
			}
			
			final Set<String> undefined_sections = new HashSet<>();
			try(BufferedReader br= super.openBufferedReader(oneFileOrNull(args))) {
				String prevComment = null;
				String line;
				Section section = null;
				while((line=br.readLine())!=null) {
					
					if (line.startsWith("# The command line was:")) {
						String[] tokens = line.split("\\s+");
						for(int i=0;i+1< tokens.length;i++) {
							if(tokens[i].equals("--depth") || tokens[i].equals("-d")) {
								tokens = CharSplitter.COMMA.split(tokens[i+1]);
								if(tokens.length==3) {
									this.depth_bin_size = Integer.parseInt(tokens[2]);
									LOG.info("depth-bin-size : "+this.depth_bin_size);
									}
								break;
								}
							}
						}
					
					if(line.startsWith("#")) {
						prevComment=line;
						continue;
						}
					final String tokens[]=CharSplitter.TAB.split(line);
					if(tokens[0].equals("ID")) {
						this.id2vcfs.put(tokens[1], String.join(",",Arrays.asList(tokens).subList(2, tokens.length)));
						continue;
						}
					if(section==null || !section.ID.equals(tokens[0])) {
						section	= sections.stream().filter(S->S.ID.equals(tokens[0])).findFirst().orElse(null);
						
						if(section==null) {
							if(undefined_sections.add(tokens[0])) {
								LOG.warn("no handler defined for "+tokens[0]);
								}
							continue;
							}
						if(!section.enabled) {
							continue;
							}
						if(prevComment==null || !section.header.equals(prevComment)) {
							LOG.error("expected "+section.header+" but got "+prevComment);
							return -1;
							}
						}
					section.add(tokens);
					}
				
				
				}
			
			try(PrintWriter pw = new PrintWriter(System.out)) {
				pw.println("# output to be piped in R");
				for(Section sec: sections) {
					if(!sec.enabled) continue;
					if(sec.lines.isEmpty()) continue;
					LOG.info("#invoking "+sec.ID);
					pw.println("# BEGIN "+sec.ID);
					sec.print(pw);
					pw.println("# END "+sec.ID);
					}
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
			new PlotBcftoolsStats().instanceMainWithExit(args);
			}
		

	}
