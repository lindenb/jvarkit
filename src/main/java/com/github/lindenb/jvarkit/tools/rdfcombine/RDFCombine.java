package com.github.lindenb.jvarkit.tools.rdfcombine;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import org.apache.jena.rdf.model.Model;
import org.apache.jena.rdf.model.ModelFactory;
import org.apache.jena.rdf.model.Statement;
import org.apache.jena.rdf.model.StmtIterator;
import org.apache.jena.vocabulary.RDF;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.IOUtil;

/**
BEGIN_DOC

## Input

input is a list of URI, or a file with the '.list' suffix containing the URI

URI are processed sequentially.

URI can be prefixed with '+' or '-' . A minus sign means that the model will be substracted to the current one.

WARNING: blank nodes are always not equal

## Example

```
$ java -jar dist/jvarkit.jar rdfcombine ${HOME}/file.rdf   -${HOME}/file.rdf 
[INFO][RDFCombine]number of statements after adding file.rdf) = 1663
[INFO][RDFCombine]number of statements after substracting -file.rdf) = 0
<rdf:RDF
    xmlns:rel="http://purl.org/vocab/relationship/"
    xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
    xmlns:bio="http://purl.org/vocab/bio/0.1/"
    xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
    xmlns:foaf="http://xmlns.com/foaf/0.1/">
</rdf:RDF>


$ java -jar dist/jvarkit.jar rdfcombine ${HOME}/file2.rdf  ${HOME}/file2.rdf 


java -jar dist/jvarkit.jar  rdfcombine "https://raw.githubusercontent.com/BruceMWhealton/Gedcom-RDF/master/Disney.rdf" "https://raw.githubusercontent.com/BruceMWhealton/Gedcom-RDF/master/Thomas.rdf"

```

END_DOC
*/
@Program(name="rdfcombine",
description="Substract/Add RDF models",
creationDate="20230903",
modificationDate="20230903",
jvarkit_amalgamion =  true
)
public class RDFCombine extends Launcher{
	
	private static final Logger LOG = Logger.of(RDFCombine.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outFile=null;

	@Parameter(names={"--base"},description="(URI) xml:base when reading rdf model from stdin or writing model")
	private String xmlBase="";
	@Parameter(names={"--keep-rdfType"},description="ignore rdf:type in minus model. (do not remove rdf:type in source)")
	private boolean keep_rdfType= false;
	@Parameter(names={"--output-format"},description="Write a serialized represention of a model in a specified language. Predefined values are RDF/XML, RDF/XML-ABBREV, N-TRIPLE, TURTLE, (and TTL) and N3")
	private String outputFormat = "RDF/XML-ABBREV";
	
	@Override
	public int doWork(final List<String> args) {
		try {
			Model sourceModel = null;
			final List<String> paths = IOUtils.unrollStrings(args);
			if(paths.isEmpty()) {
				LOG.error("empty input");
				return -1;
				}
			if(paths.stream().
				filter(F->F.equals("stdin") || F.equals("+stdin") || F.equals("-stdin")).count()>1L) {
				LOG.error("'stdin' can be used only one time");
				return -1;
				}
			
			for(int i=0;i<paths.size();i++) {
				final String uri1 = paths.get(i);
				boolean substract = false;
				final String uri2;
				if(uri1.startsWith("+")) {
					uri2 = uri1.substring(1);
					}
				else if(uri1.startsWith("-")) {
					uri2 = uri1.substring(1);
					substract = true;
					}
				else
					{
					uri2 = uri1;
					}
				if(i==0 && substract) {
					LOG.error("first rdf file "+ uri1 + " cannot be substracted");
					return -1;
					}
				final Model model  = ModelFactory.createDefaultModel();
				if(uri2.equals("stdin")) {
					if(StringUtils.isBlank(this.xmlBase)) throw new IllegalArgumentException("--base is empty");
					model.read(stdin(), this.xmlBase);
					}
				else
					{
					final String uri3;
					if(!IOUtil.isUrl(uri2)) {
						uri3 = Paths.get(uri2).toUri().toString();
						}
					else
						{
						uri3 = uri2;
						}
					model.read(uri3);
					}
				if(i==0) {
					sourceModel = model;
					}
				else if(substract) {
					final StmtIterator iter = model.listStatements();
					while(iter.hasNext()) {
						final Statement stmt = iter.next();
						if(keep_rdfType && stmt.getPredicate().equals(RDF.type)) continue;
						sourceModel.remove(stmt);
						}
					iter.close();
					}
				else
					{
					sourceModel.add(model);
					}
				LOG.info("number of statements after "+(substract?"substracting":"adding")+" "+uri1+") = "+sourceModel.size());
				}			
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outFile)) {
				if(StringUtils.isBlank(this.xmlBase)) {
					sourceModel.write(pw,this.outputFormat);
					}
				else
					{
					sourceModel.write(pw,this.outputFormat,this.xmlBase);
					}
				pw.flush();
				}
			return 0;
			}
		catch(Throwable err ) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new RDFCombine().instanceMainWithExit(args);;

	}

}
