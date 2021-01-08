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
package com.github.lindenb.jvarkit.tools.groupbygene;

import java.io.BufferedReader;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.goa.GOAFileIterator;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.go.GoTree;
import com.github.lindenb.jvarkit.util.go.GoTree.Term;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.RuntimeIOException;

public class GoGeneReporter extends Launcher{
	private static Logger LOG=Logger.build(GoGeneReporter.class).make(); 
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names="--goa",description="(goa input url)")
	private String goaUri = GOAFileIterator.DEFAULT_GOA_URI;
	@Parameter(names="--limit",description="limit output to the descendant of those GO terms. Muliple separated with spaces or comma. eg. 'GO:0015075,GO:0060047' ion transmembrane transporter activity,heart contraction")
	private String limitTermStr = "";
	@Parameter(names="--header",description="first line is header.")
	private boolean first_line_is_header= false;
	@Parameter(names="-C",description="gene column name.")
	private int geneColumnName1=1;
	@ParametersDelegate
	private com.github.lindenb.jvarkit.util.go.GoTree.ReadingGo readingGo = new  com.github.lindenb.jvarkit.util.go.GoTree.ReadingGo();


	private abstract class Reporter implements AutoCloseable{
		abstract void beginDoc();
		abstract void report(GoTree.Term term,Set<String> genes,List<List<String>> table);
		abstract void endDoc();
		}
	private class TextReporter extends Reporter{
		PrintWriter pw;
		TextReporter(PrintWriter pw) {
			this.pw=pw;
			}
		void beginDoc() {
			
		}
		void report(GoTree.Term term,Set<String> genes,List<List<String>> table) {
			pw.println(">>> " + term.getAcn()+ " : "+ term.getDefinition());
			for(int i=0;i< table.size();i++) {
				if(i==0 && first_line_is_header || genes.contains(table.get(i).get(geneColumnName1-1))) {
					pw.println(table.get(i).stream().collect(Collectors.joining("\t")));
					}
				}
			
			pw.println();
			}

		void endDoc() {
			
		}
		@Override
		public void close() throws Exception {
			pw.flush();
			pw.close();
			
		}
	}
	
	private class XmlReporter extends Reporter {
		XMLStreamWriter w;
		@Override
		void beginDoc() {
			
			}
		@Override
		void report(Term term, Set<String> genes,List<List<String>> table) {
			// TODO Auto-generated method stub
			
		}
		@Override
		void endDoc() {
			try {
				w.writeEndElement();
				w.writeEndElement();
				w.writeEndDocument();
				} 
			catch(XMLStreamException err ) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public void close() throws Exception {
			w.flush();
			w.close();
			}
		}
	
	
@Override
public int doWork(final List<String> args) {

		try {
			final String input = oneFileOrNull(args);
			final List<List<String>> table = new ArrayList<>();
			try(BufferedReader br= super.openBufferedReader(input)) {
				String line;
				while((line=br.readLine())!=null) {
					final List<String> tokens = CharSplitter.TAB.splitAsStringList(line);
					if(tokens.size()< this.geneColumnName1) {
						throw new JvarkitException.TokenErrors("expected "+this.geneColumnName1+" columns",tokens);
						}
					table.add(tokens);
					}
				}
			if(table.isEmpty()) {
				LOG.info("No data. Bye");
				return 0;
				}
			final Set<String> geneNames = table.stream().
					skip(first_line_is_header?1L:0L).
					map(T->T.get(geneColumnName1-1)).
					collect(Collectors.toSet());
			final Map<String, Set<GoTree.Term>> gene2go = new HashMap<>(geneNames.size());
			
			final GoTree mainGoTree = this.readingGo.createParser().
					setDebug(false).
					parse(this.readingGo.goUri);

			final Set<GoTree.Term> limitToTerms;
			if(StringUtils.isBlank(this.limitTermStr)) {
				limitToTerms = null;
				}	
			else
				{
				limitToTerms = Arrays.stream(this.limitTermStr.split("[ ,\t\n]+")).map(S->{
					GoTree.Term term = mainGoTree.getTermByAccession(S);
					if(term==null) term = mainGoTree.getTermByName(S);
					if(term==null) throw new IllegalArgumentException("Cannot find GO term : "+S);
					return term;
					}).collect(Collectors.toSet());
				}
			
			try(GOAFileIterator goain = GOAFileIterator.newInstance(this.goaUri)) {
				while(goain.hasNext()) {
					final GOAFileIterator.GafRecord rec = goain.next();
					if(rec.getQualifiers().contains("NOT")) continue;
					if(!geneNames.contains(rec.getObjectSymbol())) continue;
					 final GoTree.Term term =  mainGoTree.getTermByAccession(rec.getGoId());
					 if(term==null) {
						LOG.warn("Cannot find GO term "+rec.getGoId());
						continue;
					 	}
					 Set<GoTree.Term> acns = gene2go.get(rec.getObjectSymbol());
					 if(acns==null) {
						acns = new HashSet<>();
						gene2go.put(rec.getObjectSymbol(), acns);
					 	}
					 acns.add(term);
					}
				}
			LOG.warn("No GO term was found associated to the following genes:"+
				geneNames.stream().filter(G->!gene2go.containsKey(G)).collect(Collectors.joining(" ")));
			
			Reporter reporter=new TextReporter(super.openPathOrStdoutAsPrintWriter(this.outputFile));
			reporter.beginDoc();
			
			for(final GoTree.Term term: mainGoTree.getTerms()) {
				Objects.requireNonNull(term);
				if(limitToTerms!=null && limitToTerms.stream().noneMatch(T->term.isDescendantOf(T))) continue;
				final Set<String> displayGenes = gene2go.entrySet().
						stream().
						filter(KV->KV.getValue().stream().anyMatch(TERM->TERM.isDescendantOf(term))).
						map(KV->KV.getKey()).
						collect(Collectors.toSet());
				if(displayGenes.isEmpty()) continue;
				reporter.report(term,displayGenes,table);
				}
			reporter.endDoc();
			reporter.close();
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new GoGeneReporter().instanceMainWithExit(args);
	}

}
