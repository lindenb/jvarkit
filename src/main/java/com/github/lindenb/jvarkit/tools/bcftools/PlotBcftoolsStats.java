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
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

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
	/** PSC per sample count *******************************************************************/
	private class PSCSection extends Section {
		private class Column {
			final String label;
			final ToDoubleFunction<List<String>>  fun;
			Column(final String label,ToDoubleFunction<List<String>>  fun) {
				this.label= label;
				this.fun = fun;
				}
			}
		
		PSCSection() {
			super("PSC","# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions\t[8]nTransversions\t[9]nIndels\t[10]average depth\t[11]nSingletons\t[12]nHapRef\t[13]nHapAlt\t[14]nMissing");
			}
		private void print(PrintWriter w,String id,String title,String extra,List<Column> columns) {
			w.println(device(this.ID+"."+id+"."+title.toLowerCase().replace(' ', '_'),null));
			final List<List<String>> L = getLinesForId(id).stream().
					collect(Collectors.toList());
			final List<String> colNames = new ArrayList<>();
			
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
			
			w.println("T2 <- as.matrix(data.frame("+String.join(",", colNames) +"))");
				

			w.println(
					"barplot(T2,main="+quote(title)+","+
					"sub="+quote(id2vcfs.getOrDefault(id,""))+","+
	                "xlab=\"Sample\","+
	                "ylab=\"Count Variants\","+
	                "names=c("+L.stream().map(S->quote(S.get(2))).collect(Collectors.joining(","))+"),"+
	                "legend=c("+columns.stream().map(S->quote(S.label)).collect(Collectors.joining(","))+
					"),las=2"+extra+",col=rainbow("+ columns.size() +"))");
					
			w.println("dev.off()");
			}
		@Override
		void print(PrintWriter w) {
			for(String id : getIds()) {
				print(w,id,"Average Depth","",Arrays.asList(new Column("Depth",L->Double.parseDouble(L.get(9)))));
				print(w,id,"Missing","",Arrays.asList(new Column("Missing",L->Double.parseDouble(L.get(13)))));
				print(w,id,"Genotypes","",
						Arrays.asList(
								new Column("HOM_REF",L->Double.parseDouble(L.get(3))),
								new Column("HET",L->Double.parseDouble(L.get(5))),
								new Column("HOM_VAR",L->Double.parseDouble(L.get(4))),
								new Column("NO_CALL",L->Double.parseDouble(L.get(13)))
								));
				}
			}
		}
	/* https://stackoverflow.com/a/52498075 */
	private String getColorsForSamples(final List<String> samples) {
		if(this.phenotype2samples.size()<2) return "";
		final String hue = samples.stream().
			map(S->this.sample2phenotype.get(S)).
			map(P->{
				double n= this.phenotype2samples.size();
				int i=0;
				for(String k:this.phenotype2samples.keySet()) {
					if(k.equals(P)) break;
					i++;
					}
				return String.valueOf(i/n);
				}
			).collect(Collectors.joining(","));
		return "hsv(h=c("+hue+"))";
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
				new DPSection(),
				new STSection(),
				new QualSection(),
				new SNSection(),
				new PSCSection(),
				new DPSection()
				);
			
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
