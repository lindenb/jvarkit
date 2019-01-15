/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.pubmed;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;

/**
BEGIN_DOC

## Example

```
$ java -jar dist/pubmeddump.jar 'Bioinformatics 2017' | java -jar dist/pubmedcodinglang.jar

(...)
28464826	java	microTaboo: a general and practical solution to the k-disjoint problem.	2017	, running [java] 7 and hig
28453684	python	MAPPI-DAT: data management and analysis for protein-protein interaction data from the high-throughput MAPPIT cell microarray platform.	2017	eloped in [python], using r 
28453672	python	Primerize-2D: automated primer design for RNA multidimensional chemical mapping.	2017	and-alone [python] package f
28449031	python	RNAblueprint: Flexible multiple target nucleic acid sequence design.	2017	ritten in [python] to demons
28444126	ruby	CImbinator: A web-based tool for drug synergy analysis in small- and large-scale datasets.	2017	ritten in [ruby] and r. it
28437135	python	D3GB: An Interactive Genome Browser for R, Python, and WordPress.	2017	ackage, a [python] module, a
28431529	python	HirBin: high-resolution identification of differentially abundant functions in metagenomes.	2017	nted as a [python] package a
28430977	java	Deep Mining Heterogeneous Networks of Biomedical Linked Data to Predict Novel Drug-Target Associations.	2017	eloped in [java] and it is
28415074	python	NaviCom: a web application to create interactive molecular network portraits using multi-level omics data.	2017	avicom, a [python] package a

```

Plotting:

```
cut -f2,4 output.txt  | sort | uniq -c |\
awk '{printf("%s\t%s\t%s\n",$3,$2,$1);}' | ./yxv2table -n 0 | grep -v '^null' | grep -v '^     ' | sed 's/^#/Year/' > table.txt && \
gnuplot script.gnuplot
``` 



```
set terminal png size 2000,1500;
set output "output.png"
set title "Bioinformatics/language"
set key invert reverse Left outside
set key autotitle columnheader
set yrange [0:]
set auto x
unset xtics
set xlabel "Year"
set ylabel "Count"
set xtics nomirror rotate by -45 scale 0
set style data histogram
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.75
N=`awk 'NR==1 {print NF}' table.txt`
plot 'table.txt' using 2:xtic(1), for [i=3:N] '' using i;
```


![https://pbs.twimg.com/media/C_AZTzFXoAUpiuI.jpg:large](https://pbs.twimg.com/media/C_AZTzFXoAUpiuI.jpg:large)

## See also

* https://gist.github.com/lindenb/83196adbb034ef5874086d10dd9772ac 

END_DOC
 */
@Program(name="pubmedcodinglang",
	description="Programming language use distribution from recent programs / articles",
	keywords={"pubmed","xml","code","programming"},
	biostars=251002
	)
public class PubmedCodingLanguages
	extends Launcher
	{
	private static final Logger LOG = Logger.build(PubmedCodingLanguages.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outFile=null;
	@Parameter(names={"-c","--common"},description="[20180302] only common languages. What is a 'common' language ? well .. it's subjective...")
	private boolean common_only_false;

	private interface ProgLanguage {
		/** language name */
		public String getName();
		/** context */
		public String match(final String s);
		// is a commong language ?
		public boolean isCommon();
	}
	
	private static class ProgLanguageImpl implements ProgLanguage
		{
		protected final Set<String> othersNames;
		private final String label;
		private boolean commonFlag = true;
		ProgLanguageImpl( final String name,final String...others) {
			this.label = name.toLowerCase();
			this.othersNames = new LinkedHashSet<>();
			this.othersNames.add(name);
			this.othersNames.addAll(Arrays.asList(others));
		}
		@Override
		public String getName() {
			return label;
		}
		private boolean isDelim(char c) {
			return c==',' || c=='.' || Character.isWhitespace(c) || c=='\'';
		}
		@Override public String match(final String s) {
			return this.othersNames.
					stream().
					map(L->search(L,s)).
					filter(S->!StringUtil.isBlank(S)).
					findFirst().orElse(null);
					
		}
		ProgLanguageImpl flagHasNotCommon() {
			this.commonFlag= false;
			return this;
		}
		@Override
		public boolean isCommon() {
			return this.commonFlag;
			}
		
		private String search(final String key,final String s) {
			int i=0;
			if(s==null || s.isEmpty()) return null;
			while((i=s.indexOf(key.toLowerCase(),i))!=-1) {
				int j = i+key.length();
				boolean test_end = (j>=s.length() || isDelim(s.charAt(j)));
				boolean test_start = (i==0 || isDelim(s.charAt(i-1)));
				if( test_end && test_start)
					{
					final int extend=10;
					return s.substring(Math.max(0,i-extend),i)+
							"["+key+"]"+
							s.substring(j,Math.min(s.length(),j+extend));
					}
				++i;
				}

			return null;
			}
		}
	
	
	
	private final List<ProgLanguage> languages = new ArrayList<>();
	
		
	public PubmedCodingLanguages()
		{
		this.languages.add(new ProgLanguageImpl("java","biojava","jre","jdk","jsp","javac","j2ee","webstart") {
			@Override
			public String match(String s) {
				s = s.replaceAll(".jsp", " ");
				return super.match(s);
			}
		});
		
		this.languages.add(new ProgLanguageImpl("ruby","rails").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("python","cython","django","python3","python2.7","biopython"));
		this.languages.add(new ProgLanguageImpl("awk","nawk","gawk").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("golang","go language","go programming").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("couchdb","couch database").flagHasNotCommon());
		//this.languages.add(new ProgLanguageImpl("R"));
				
		this.languages.add(new ProgLanguageImpl("C++") {
			@Override
			public String match(String s) {
				String m = super.search("C++", s);if(m!=null) return m;
				m = super.search("C ++", s);if(m!=null) return m;
				m = super.search("c++14", s);if(m!=null) return m;
				return null;
			}
		});

		
		this.languages.add(new ProgLanguageImpl("C","lex/yacc","flex/bison") {
			@Override
			public String match(String s) {
				s=s.toLowerCase().replaceAll("c. elegans", " ");
				String m = super.search("written in C", s);if(m!=null) return m;
				m = super.search("C language", s);if(m!=null) return m;
				m = super.search("language C", s);if(m!=null) return m;
				m = super.search("C programming", s);if(m!=null) return m;
				m = super.search("'C'", s);if(m!=null) return m;
				m = super.search("using C", s);if(m!=null) return m;
				m = super.search("ansi-C", s);if(m!=null) return m;
				return null;
			}
		});
		
		this.languages.add(new ProgLanguageImpl("R","r-shiny") {
			@Override
			public String match(String s) {
				String m = super.search("R package", s);if(m!=null) return m;
				m = super.search("R language", s);if(m!=null) return m;
				m = super.search("language r", s);if(m!=null) return m;
				
				m = super.search("Rstat", s);if(m!=null) return m;
				m = super.search("'R'", s);if(m!=null) return m;
				m = super.search("using R", s);if(m!=null) return m;
				m = super.search("written in R", s);if(m!=null) return m;
				m = super.search("and R", s);if(m!=null) return m;
				m = super.search("CRAN", s);if(m!=null) return m;
				m = super.search("Bioconductor", s);if(m!=null) return m;
				return null;
			}
		});

		this.languages.add(new ProgLanguageImpl("scala").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("makefile","gnu make").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("mongodb").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("sql").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("sqlite","sqlite3").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("mysql").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("bash","shell script"));
		this.languages.add(new ProgLanguageImpl("perl","CPAN","bioperl"));

		
		this.languages.add(new ProgLanguageImpl("php") {
			@Override
			public String match(String s) {
				s = s.replaceAll(".php", " ");
				return super.match(s);
				}
			});
		this.languages.add(new ProgLanguageImpl("c#").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("objective-c").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("haskell").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("lua").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("clojure").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("mathlab","matlab"));
		this.languages.add(new ProgLanguageImpl("groovy").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("javascript"));
		this.languages.add(new ProgLanguageImpl("nodejs").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("pascal").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("erlang").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("smalltalk").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("fortran","f77"));
		this.languages.add(new ProgLanguageImpl("macro excel","excel macros","excel macro","spreadsheet excel","microsoft excel").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("microsoft access").flagHasNotCommon());
		this.languages.add(new ProgLanguageImpl("visual basic").flagHasNotCommon());
		
		}
	
	/** extract text from stream. Cannot use XMLEventReader.getTextContent() 
	 * when a public title contains some tag like '<sup>'
	 */
	private String textContent(final XMLEventReader r) throws XMLStreamException
		{
		final StringBuilder sb=new StringBuilder();
		while(r.hasNext())
			{
			final XMLEvent evt = r.nextEvent();
			if(evt.isEndElement()) break;
			else if(evt.isStartElement())
				{
				sb.append(textContent(r));
				}
			else if(evt.isCharacters())
				{
				sb.append(evt.asCharacters().getData());
				}
			}
		return sb.toString();
		}
	
	private void scanArticle(
			PrintStream out,
			final String rootName,
			final XMLEventReader r
			) throws XMLStreamException,IOException {
			String article_pmid=null;
			String article_title="" ;
			String article_year=null;
			String abstractText = "";
			boolean PubDate=false;
			while(r.hasNext()) {
				final XMLEvent evt=r.nextEvent();
				
				if(evt.isStartElement()) {
					final StartElement start = evt.asStartElement();
					final String eltName = start.getName().getLocalPart();
					if(article_pmid==null && eltName.equals("PMID")) {
						article_pmid=r.getElementText();
					}
					else if(article_title.isEmpty() && eltName.equals("ArticleTitle")) {
						article_title= textContent(r);
					}
					else if(eltName.equals("PubDate")) {
						PubDate=true;
					}
					else if(article_year==null && PubDate && eltName.equals("Year")) {
						article_year=r.getElementText();
					}
					else if(eltName.equals("Abstract")) {
						abstractText = textContent(r);
						}
					}
			else if(evt.isEndElement()) {
				final EndElement end = evt.asEndElement();
				final String eltName = end.getName().getLocalPart();
				if(eltName.equals("PubDate")) {
					PubDate=false;
					}
				else if(eltName.equals(rootName)) 
					{
					break;
					}
				}
			
			}//end of xml read
		final String normalizedAbstract = (article_title+" "+abstractText).replace('\n', ' ').toLowerCase();
		boolean found_something = false;
		for(final ProgLanguage pg: this.languages)
			{
			if(this.common_only_false && !pg.isCommon()) continue;
			final String context = pg.match(normalizedAbstract);
			if(context==null) continue;
			out.print(article_pmid);
			out.print('\t');
			out.print(pg.getName());
			out.print('\t');
			out.print(article_title);
			out.print('\t');
			out.print(article_year);
			out.print('\t');
			out.print(context);
			out.println();
			found_something = true;
			}
		if(!found_something && 
			(normalizedAbstract.contains("programming language") || normalizedAbstract.contains("written in")))
			{
			LOG.warning("Cannot get information in:\t"+normalizedAbstract);
			}
		}	
	@Override
	public int doWork(final List<String> args) {
		final String inputName= oneFileOrNull(args);
		PrintStream out=null;
		XMLEventReader r=null;
		InputStream in=null;
		try {
			final XMLInputFactory xmlInputFactory = XMLInputFactory.newFactory();
			xmlInputFactory.setXMLResolver(new XMLResolver() {
				@Override
				public Object resolveEntity(String publicID, String systemID, String baseURI, String namespace)
						throws XMLStreamException {
					LOG.debug("Ignoring resolve Entity");
					return new ByteArrayInputStream(new byte[0]);
				}
			});
			in=(inputName==null?stdin():IOUtils.openURIForReading(inputName));
			r = xmlInputFactory.createXMLEventReader(in);
			
			out = super.openFileOrStdoutAsPrintStream(this.outFile);
			while(r.hasNext()) {
				final XMLEvent evt= r.nextEvent();
				if(evt.isStartElement() )	{
					final String localName= evt.asStartElement().getName().getLocalPart();
					if(localName.equals("PubmedArticle") || localName.equals("PubmedBookArticle"))
						{
						scanArticle(out,localName,r);	
						}
					}
				}
			
			r.close();r=null;
			in.close();in=null;
			out.flush();out.close();out=null;
			return 0;
		} catch (final Exception err) {
			LOG.error(err);
			return -1;
		} finally {
			CloserUtil.close(r);
			CloserUtil.close(in);
			CloserUtil.close(out);
		}

		}
		
	
	public static void main(String[] args) {
		new PubmedCodingLanguages().instanceMainWithExit(args);
	}
}
