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

*/
package com.github.lindenb.jvarkit.tools.vcfsparql;

import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.apache.jena.query.Query;
import org.apache.jena.query.QueryExecution;
import org.apache.jena.query.QueryExecutionFactory;
import org.apache.jena.query.QueryFactory;
import org.apache.jena.query.ResultSet;
import org.apache.jena.query.ResultSetFormatter;
import org.apache.jena.rdf.model.Model;
import org.apache.jena.rdf.model.ModelFactory;
import org.apache.jena.sparql.resultset.ResultsFormat;
import org.apache.jena.sys.JenaSubsystemRegistry;
import org.apache.jena.sys.JenaSubsystemRegistryBasic;
import org.apache.jena.sys.JenaSystem;
import org.apache.jena.vocabulary.DC;
import org.apache.jena.vocabulary.RDF;
import org.apache.jena.vocabulary.RDFS;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

/**
 
BEGIN_DOC
 
## Example
 
``` 
$ cat query.sparql


PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>


SELECT distinct *
WHERE {
    ?s ?p ?o .
    }
```

```
$ java -jar dist/vcfsparql.jar -f query.sparql src/test/resources/rotavirus_rf.vcf.gz

---------------------------------------------------------------------------------------------------------------
| s      | p                                                 | o                                              |
===============================================================================================================
| _:b0   | <jvarkit:allele>                                  | "C"                                            |
| _:b0   | <jvarkit:sample>                                  | "S5"                                           |
| _:b1   | <jvarkit:type>                                    | "HOM_VAR"                                      |
| _:b1   | <jvarkit:genotype>                                | _:b0                                           |
| _:b0   | <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> | <jvarkit:Genotype>                             |
| _:b2   | <jvarkit:allele>                                  | "A"                                            |
| _:b2   | <jvarkit:sample>                                  | "S4"                                           |
| _:b1   | <jvarkit:type>                                    | "HOM_REF"                                      |
| _:b1   | <jvarkit:genotype>                                | _:b2                                           |
| _:b2   | <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> | <jvarkit:Genotype>                             |
| _:b3   | <jvarkit:allele>                                  | "A"                                            |
| _:b3   | <jvarkit:sample>                                  | "S3"                                           |
| _:b1   | <jvarkit:genotype>                                | _:b3                                           |
| _:b3   | <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> | <jvarkit:Genotype>                             |
| _:b4   | <jvarkit:allele>                                  | "A"                                            |
| _:b4   | <jvarkit:sample>                                  | "S2"                                           |
| _:b1   | <jvarkit:genotype>                                | _:b4                                           |
| _:b4   | <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> | <jvarkit:Genotype>                             |
| _:b5   | <jvarkit:allele>                                  | "A"                                            |
| _:b5   | <jvarkit:sample>                                  | "S1"                                           |
(...)
---------------------------------------------------------------------------------------------------------------
```

 
END_DOC
 
 
 */

@Program(name="vcfsparql",
description="Query RDF with Sparql. Very slow. Just a proof of concept",
keywords={"vcf","sparql","rdf","arq","semanticweb"},
creationDate="2019-03-06",
modificationDate="2019-03-06"
)
public class VcfSparql extends Launcher {
	private static final Logger LOG=Logger.build(VcfSparql.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-f","--sparql-file"},description="Query as sparql file")
	private Path sparqlQueryFile = null;
	@Parameter(names={"-e","--sparql"},description="Query as sparql string")
	private String sparqlQueryString  = null;
	@Parameter(names={"-a","--prefixes"},description="Prepend common RDF prefixes to the expression")
	private boolean add_prefixes = false;
	@Parameter(names={"--format"},description="output format. One of https://jena.apache.org/documentation/javadoc/arq/org/apache/jena/sparql/resultset/ResultsFormat.html#lookup-java.lang.String-")
	private String outputFormat = "text";
	@Parameter(names={"--code"},description="show code")
	private boolean showCode = false;
	@Parameter(names={"-r","--region"},description="limit query to this genomic interval. "+IntervalParser.OPT_DESC)
	private String regionStr = null;

	@Override
	public int doWork(final List<String> args) {
		List<VariantGraph> graphs = new ArrayList<>();
		Model model = null;
		OutputStream out = stdout();
		if(sparqlQueryFile==null && StringUtils.isBlank(this.sparqlQueryString))
			{
			LOG.error("SPARQL query undefined");
			return -1;
			}
		else if(sparqlQueryFile!=null && !StringUtils.isBlank(this.sparqlQueryString))
			{
			LOG.error("SPARQL query expression and file both undefined");
			return -1;
			}
	
		
		try {
			// https://jena.apache.org/documentation/notes/system-initialization.html
			JenaSubsystemRegistry r = new JenaSubsystemRegistryBasic() {
			    @Override
			    public void load() {
			        if ( JenaSystem.DEBUG_INIT ) LOG.debug("Custom load") ;
			        this.add(new org.apache.jena.riot.system.InitRIOT());
			        this.add(new org.apache.jena.sparql.system.InitARQ());
			    }
			} ;

			// Set the sub-system registry
			JenaSystem.setSubsystemRegistry(r);
			JenaSystem.init();
			
			final  ResultsFormat rFmt = ResultsFormat.lookup(this.outputFormat);
			if(rFmt==null) {
				LOG.info("bad output format : "+ this.outputFormat);
				return -1;
				}

			
			IOUtils.unrollFiles2018(args);
			
			String queryString = this.sparqlQueryFile!=null ?
					new String(Files.readAllBytes(this.sparqlQueryFile)):
					this.sparqlQueryString
					;
			if(StringUtils.isBlank(queryString)) {
				LOG.error("query is empty");
				return -1;
				}
			if(add_prefixes) {
				queryString = 
					"PREFIX rdf: <" +RDF.getURI()+">\n" +
					"PREFIX rdfs: <" +RDFS.getURI()+">\n" +
					"PREFIX dc: <" +DC.getURI()+">\n" +
					"PREFIX vcf: <" +VariantGraph.NS+">\n" +
					"\n" + queryString;
				}
			
			for(final Path vcfInput : IOUtils.unrollPaths(args)) {
				final VariantGraph graph = new VariantGraph(vcfInput);
				
				if(!StringUtils.isBlank(this.regionStr)) {
					final IntervalParser parser= new IntervalParser(SAMSequenceDictionaryExtractor.extractDictionary(vcfInput));
					final Interval rgn = parser.parse(this.regionStr);
					if(rgn==null) {
						LOG.error("cannot parse interval "+this.regionStr);
						return -1;
						}
					graph.setInterval(rgn);
					}
				
				graphs.add(graph);
				}
			
			if(graphs.isEmpty())
				{
				LOG.error("no vcf file defined");
				return -1;
				}	
			
			for(int i=0;i< graphs.size();i++) {
				final Model tmpModel = ModelFactory.createModelForGraph(graphs.get(i));
				if(i==0) {
					model = tmpModel;
					}
				else
					{
					model = ModelFactory.createUnion(model, tmpModel);
					}
				}
			
			if(this.showCode) {
				LOG.info("\n"+queryString+"\n");
			}
			final Query query= QueryFactory.create(queryString,"zzz");
			final QueryExecution execution = QueryExecutionFactory.create(query, model);
			ResultSet results =execution.execSelect();
			
			
			out= super.openPathOrStdoutAsStream(this.outputFile);
			ResultSetFormatter.output(out,results, rFmt);
			out.flush();
			out.close();
			out=null;
			
			model.close();
			model=null;
			graphs.stream().forEach(G->G.close());graphs.clear();
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(model);
			graphs.stream().forEach(G->G.close());
			CloserUtil.close(out);
			}
		
		}
	
	public static void main(String[] args) {
		new VcfSparql().instanceMainWithExit(args);
	}
}
