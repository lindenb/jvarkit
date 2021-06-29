package com.github.lindenb.jvarkit.tools.bcftools;

import java.io.BufferedReader;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.RangeOfDoubles;
import com.github.lindenb.jvarkit.util.Counter;
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
	@DynamicParameter(names = "-D", description = "set some css style elements. '-Dkey=value'. Undocumented.")
	private Map<String, String> dynaParams = new HashMap<>();


	private enum Format {PDF,PNG,SVG};
	private final Map<String,String> id2file = new HashMap<>();

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
					"sub="+quote(id2file.getOrDefault(id,""))+","+
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
					"sub="+quote(id2file.getOrDefault(id,""))+","+
	                "ylab=\"Count Variants\","+
	                "names=c("+L.stream().map(T->quote(removeColumn(T.get(2)))).collect(Collectors.joining(","))+"),"+
					"las=2,cex.names=0.5)");
					}

				w.println("dev.off()");
				}
			}
	/** PSC *******************************************************************/
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
					"sub="+quote(id2file.getOrDefault(id,""))+","+
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
				new QualSection(),
				new SNSection(),
				new PSCSection()
				);
			final Set<String> undefined_sections = new HashSet<>();
			try(BufferedReader br= super.openBufferedReader(oneFileOrNull(args))) {
				String prevComment = null;
				String line;
				Section section = null;
				while((line=br.readLine())!=null) {
					if(line.startsWith("#")) {
						prevComment=line;
						continue;
						}
					final String tokens[]=CharSplitter.TAB.split(line);
					if(tokens[0].equals("ID")) {
						id2file.put(tokens[1], tokens[2]);
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
