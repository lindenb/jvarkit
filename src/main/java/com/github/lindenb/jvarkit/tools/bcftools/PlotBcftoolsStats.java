package com.github.lindenb.jvarkit.tools.bcftools;

import java.io.BufferedReader;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
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


	private enum Format {PDF,PNG,SVG};

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
	
	private class QualSection extends Section {
		QualSection() {
			super("QUAL","# QUAL\t[2]id\t[3]Quality\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)\t[7]number of indels");
			}
		@Override
			void print(PrintWriter w) {
				for(String id : getIds()) {
					w.println(device(this.ID+"."+id,null));
					List<List<String>> L = getLinesForId(id);
					final double median = 1000;
					/*2.0 * L.stream().
						mapToDouble(T->Long.parseLong(T.get(3))*Double.parseDouble(T.get(2))).
						average().
						orElse(0.5);*/
					L = L.stream().
						filter(T->Double.parseDouble(T.get(2))<=median).
						collect(Collectors.toList());

					final long max_num = L.stream().mapToLong(T->Long.parseLong(T.get(3))).max().orElse(1L);
					final String xylim= "xlim=c("+
						L.stream().mapToDouble(T->Double.parseDouble(T.get(2))).min().orElse(0.0)+","+
						L.stream().mapToDouble(T->Double.parseDouble(T.get(2))).max().orElse(0.0)+
						"),ylim=c(0,"+max_num+")";
						;
					w.println("x_data<-c("+L.stream().map(T->T.get(2)).collect(Collectors.joining(","))+")");
					w.println("colors <- rainbow(4)");
					for(int i=0;i< 4;++i) {
							final int final_i = i;
							w.println("y_data<-c("+L.stream().map(T->T.get(3+final_i)).collect(Collectors.joining(","))+")");
							w.println("T2<-as.matrix(data.frame(x_data,y_data))");

							if(i==0) {
							w.println(
								"plot(T2,type = \"p\","
								+ "main=\"Stats by quality\","+
				                "xlab=\"Quality\","+
				                "ylab=\"Count\","+
				                xylim+","+
				                "col= colors[1],"+
				                "pch=3"+
				                ")");
							}
						else {
							w.println("par(new=TRUE)");
							w.println("plot(T2,type = \"p\",axes=FALSE,ann=FALSE," +
								"col= colors["+(i+1)+"],"+ xylim+",pch=3)");
							}
						}
					w.println("legend(\"topright\",legend=c("+
							Arrays.stream(CharSplitter.TAB.split(this.header)).skip(3).map(S->quote(removeColumn(S))).collect(Collectors.joining(","))+
							"),pch=16,col=colors)");

					w.println("dev.off()");
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
				new QualSection()
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
					sec.print(pw);
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
