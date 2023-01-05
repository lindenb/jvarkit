/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.hpo;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.obo.OBOParser;
import com.github.lindenb.jvarkit.obo.OBOntology;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Codec.DecodeDepth;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.gff.Gff3Writer;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;

/**
BEGIN_DOC

## Example

```
$ java -jar dist/hpoutils.jar -A "HP:0001626" | head
#ACN	NAME	DEFINITION
HP:0001195	Single umbilical artery	Single umbilical artery (SUA) is the absence of one of the two umbilical arteries surrounding the fetal bladder and in the fetal umbilical cord.
HP:0001136	Retinal arteriolar tortuosity	The presence of an increased number of twists and turns of the retinal arterioles.
HP:0025188	Retinal vasculitis	Inflammation of retinal blood vessels as manifested by perivascular sheathing or cuffing, vascular leakage and/or occlusion.
HP:0025169	Left ventricular systolic dysfunction	Abnormality of left ventricular contraction, often defined operationally as an ejection fraction of less than 40 percent.
HP:0025168	Left ventricular diastolic dysfunction	Abnormal function of the left ventricule during left ventricular relaxation and filling.
HP:0410173	Increased circulating troponin I concentration	An increased concentration of tropnin I in the blood, which is a cardiac regulatory protein that controls the calcium mediated interaction between actin and myosin. Raised cardiac troponin concentrations are now accepted as the standard biochemical marker for the diagnosis of myocardial infarction.
HP:0410174	Increased circulating troponin T concentration	An increased concentration of tropnin T in the blood, which is a cardiac regulatory protein that controls the calcium mediated interaction between actin and myosin. Raised cardiac troponin concentrations are now accepted as the standard biochemical marker for the diagnosis of myocardial infarction.
HP:0410267	Intestinal hemangioma	A hemangioma, a benign tumor of the vascular endothelial cells, located in the intestines, which includes the bowel.
HP:0410268	Spleen hemangioma	A hemangioma, a benign tumor of the vascular endothelial cells, that is located in the spleen.
```

## Example

````
$ java -jar dist/hpoutils.jar -A "HP:0001626" -action genes  -g2p genes_to_phenotype.txt  |head
8195	MKKS	HP:0000822	Hypertension	-	HP:0040282		orphadata	ORPHA:110
8195	MKKS	HP:0001629	Ventricular septal defect	-	HP:0040283		orphadata	ORPHA:2473
8195	MKKS	HP:0001631	Atrial septal defect	-	HP:0040283		orphadata	ORPHA:2473
8195	MKKS	HP:0004383	Hypoplastic left heart	-	HP:0040283		orphadata	ORPHA:2473
8195	MKKS	HP:0001636	Tetralogy of Fallot	-	HP:0040283		orphadata	ORPHA:2473
8195	MKKS	HP:0001643	Patent ductus arteriosus	-	HP:0040283		orphadata	ORPHA:2473
8195	MKKS	HP:0030680	Abnormality of cardiovascular system morphology	-		-	mim2gene	OMIM:236700
90121	TSR2	HP:0001629	Ventricular septal defect	-	HP:0040283		orphadata	ORPHA:124
90121	TSR2	HP:0001631	Atrial septal defect	-	HP:0040283		orphadata	ORPHA:124
90121	TSR2	HP:0001680	Coarctation of aorta	-	HP:0040284		orphadata	ORPHA:124
```

END_DOC
 *
 */
@Program(
		name="hpoutils",
		description="Human Phenotype Ontology Utils.",
		keywords={"phenotype","ontology","hpo","hpoa"},
		creationDate="20230105",
		modificationDate="20230105"
		)
public class HpoUtils
	extends Launcher
	{
	private static Logger LOG=Logger.build(HpoUtils.class).make();
	private enum Action{dump_table,genes,gff3};
		
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-action","--action"},description= "What shoud I do ? default is dump as table.")
	private Action action = Action.dump_table;
	@Parameter(names={"-hpo","--hpo","--obo"},description="HPO ontology in OBO format")
	private String hpoURI = "https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo";
	@Parameter(names={"-g2p","--g2p"},description="HPO gene to phenotype file. eg://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt ( Format: entrez-gene-id<tab>entrez-gene-symbol<tab>HPO-Term-ID<tab> )")
	private Path gene_to_phenotype_uri = null;
	@Parameter(names={"-g","--gff","--gff3"},description="GFF3 file for action=genes or action=gff3.")
	private String gffPath = null;
	
	@Parameter(names= {"-A","--accession",},description="User Go Terms accession numbers or name.eg HP:0001626 Abnormality of the cardiovascular system ")
	private Set<String> userAccStrings = new HashSet<>();
	@Parameter(names= {"--exclude-accession","--exclude"},description="User Go Terms to be EXCLUDED accession numbers or name.eg")
	private Set<String> excludeUserAccStrings = new HashSet<>();
	@Parameter(names= {"-af","--accession-file"},description="File containing accession numbers. One per line.")
	private Path accessionFile=null;
	@Parameter(names= {"-i","--inverse",},description="inverse the result")
	private boolean inverse=false;
	@Parameter(names= {"--debug",},description="debug",hidden=true)
	private boolean do_debug =false;
    
	private OBOntology mainHpoTree=null;
	
	
	private static class GeneToPheno {
		final String[] tokens;
		GeneToPheno(final String[] tokens) {
			this.tokens = tokens;
		}
		
		private String get(int i) {
			if(i<0 || i>=tokens.length ) {
				throw new ArrayIndexOutOfBoundsException("i="+i+" "+this.toString());
			}
			return this.tokens[i];
		}
		
		public String getGene() {
			return get(1);
			}
		public String getHpoAcn() {
			return get(2);
			}
		@Override
		public String toString() {
			return String.join("\t",this.tokens);
			}
		}
	
	private static class GeneToPhenoIterator extends AbstractCloseableIterator<GeneToPheno> {
		private final BufferedReader br;
		GeneToPhenoIterator(final Path uri) throws IOException {
			this.br = IOUtils.openPathForBufferedReading(uri);
			}
		@Override
		protected GeneToPheno advance() {
			String line;
			try {
				for(;;) {
					line = br.readLine();
					if(line==null) return null;
					if(line.startsWith("#")) continue;
					return new GeneToPheno(CharSplitter.TAB.split(line));
					}
			} catch (IOException e) {
				throw new RuntimeIOException(e);
				}
			}
		@Override
		public void close() {
			try {this.br.close();} catch(IOException err) {}
			}
		}
	
	
	public HpoUtils() {
	}
	
	
	
	private void dumpGff3(Gff3Writer out,Gff3Feature feat,final Set<String> geneNames) throws IOException {
		if(feat==null) return;
		if(feat.getType().equals("gene")) {
			if(feat.getAttribute("Name").stream().anyMatch(S->geneNames.contains(S))) {
				out.addFeature(feat);
				for(final Gff3Feature c:feat.getDescendents())  {
					out.addFeature(c);
					}
				}
			}
		else if(feat.hasChildren()) {
			for(Gff3Feature c:feat.getChildren())  {
				dumpGff3(out, c, geneNames);
				}
		}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			this.mainHpoTree = new OBOParser().
					parseObo(this.hpoURI);
			
			final Set<OBOntology.Term> userTerms = new HashSet<>();
			final Set<OBOntology.Term> excludeUserTerms = new HashSet<>();
			for(final String s:this.userAccStrings.stream().flatMap(S->Arrays.stream(S.split("[ ,;]"))).collect(Collectors.toSet()))
				{
				if(StringUtil.isBlank(s)) continue;
				final OBOntology.Term t= this.mainHpoTree.findTermByAcnOrName(s);
				if(t==null)
					{
					LOG.error("cannot find user term \""+s+"\"");
					return -1;
					}
				userTerms.add(t);
				}
			for(final String s:this.excludeUserAccStrings.stream().flatMap(S->Arrays.stream(S.split("[ ,;]"))).collect(Collectors.toSet()))
				{
				if(StringUtil.isBlank(s)) continue;
				final OBOntology.Term t= this.mainHpoTree.findTermByAcnOrName(s);
				if(t==null)
					{
					LOG.error("cannot find user excluded term \""+s+"\"");
					return -1;
					}
				excludeUserTerms.add(t);
				}
			
			
			final Predicate<OBOntology.Term> keepTerm = T->{
				boolean keep=false;
				if(userTerms.isEmpty())
					{
					keep=true;
					}
				else if(userTerms.
						stream().
						anyMatch(USERTERM->T.isDescendantOf(USERTERM))) {
					keep=true;
					}
				if(keep && !excludeUserTerms.isEmpty()) {
					 if(excludeUserTerms.
								stream().
								anyMatch(USERTERM->T.isDescendantOf(USERTERM))) {
						keep = false;
					 	}
					}
				
				if(this.inverse) keep=!keep;
				return keep;
				};
			
			if(this.accessionFile!=null)
				{
				try(BufferedReader r=IOUtils.openPathForBufferedReading(this.accessionFile)) {
					String line;
					while((line=r.readLine())!=null) {
						if(line.isEmpty() || line.startsWith("#")) continue;
						int last=0;
						for(last=0;last< line.length();++last) {
							if(Character.isWhitespace(line.charAt(last))) break;
							}
						final String s= line.substring(0, last);
						OBOntology.Term t=  this.mainHpoTree.findTermByAcnOrName(s);
						if(t==null)
								{
								LOG.error("In "+this.accessionFile+" cannot find user term \""+s+"\"");
								return -1;
								}
						userTerms.add(t);
						}
					}
				}
			
			
			switch(this.action)
				{
				case genes:
					{
					if(!args.isEmpty())
						{
						LOG.error("too many arguments");
						return -1;
						}
					if(this.gene_to_phenotype_uri==null) {
						LOG.error("undefined gene_to_phenotype_uri");
						return -1;
						}
					
					final Set<String> acns_set = this.mainHpoTree.
						stream().
						filter(keepTerm).
						map(T->T.getAcn()).
						collect(Collectors.toSet());
						
					try(GeneToPhenoIterator goain= new GeneToPhenoIterator(this.gene_to_phenotype_uri)) {
						try(PrintWriter out = super.openFileOrStdoutAsPrintWriter(this.outputFile)) {
							while(goain.hasNext()) {
								final GeneToPheno rec = goain.next();
								if(!acns_set.contains(rec.getHpoAcn())) continue;
								out.println(rec.toString());
								}
							out.flush();
							}
						}
				
					break;
					}
				case gff3:
					{
					if(this.gene_to_phenotype_uri==null) {
						LOG.error("undefined gene to phenotype uri");
						return -1;
						}
					final String input;
					if(!StringUtils.isBlank(this.gffPath)) {
						input = oneFileOrNull(args);
						}
					else {
						input = this.gffPath;
						}
					final Set<String> acns_set = this.mainHpoTree.
						stream().
						filter(keepTerm).
						map(T->T.getAcn()).
						collect(Collectors.toSet());
					final Set<String> geneNames = new HashSet<>();
					try(GeneToPhenoIterator goain= new GeneToPhenoIterator(this.gene_to_phenotype_uri)) {
							while(goain.hasNext()) {
								final GeneToPheno rec = goain.next();
								if(!acns_set.contains(rec.getHpoAcn())) continue;
								geneNames.add(rec.getGene());
								}
							
						}
					
					final Gff3Codec gff3 = new Gff3Codec(DecodeDepth.DEEP);
					try(InputStream is =(input == null? stdin():IOUtils.openURIForReading(input))) {
						final AsciiLineReader asciiLineReader = AsciiLineReader.from(is);
						final LineIterator lr= new LineIteratorImpl(asciiLineReader);
						try(OutputStream out = super.openFileOrStdoutAsStream(this.outputFile)) {
							Gff3Writer gw = new Gff3Writer(out);
							while(!gff3.isDone(lr)) {
								dumpGff3(gw,gff3.decode(lr),geneNames);
								}
							out.flush();
							}
						gff3.close(lr);
						asciiLineReader.close();
						}
					break;
					}
				case dump_table://through
				default:
					{
					if(!args.isEmpty()) {
						LOG.error("too many arguments");
						return -1;
						}
					try(PrintWriter out = super.openFileOrStdoutAsPrintWriter(this.outputFile)) {
						out.println("#ACN\tNAME\tDEFINITION");
						for(final OBOntology.Term t:this.mainHpoTree)
							{
							if(keepTerm.test(t))
								{
								out.print(t.getAcn());
								out.print('\t');
								out.print(t.getLabel());
								out.print('\t');
								out.print(t.getDescription().orElse("."));
								out.println();
								}
							}
						out.flush();
						}
					break;
					}
				}
				return 0;
				}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}

	public static void main(final String[] args)
		{
		new HpoUtils().instanceMainWithExit(args);
		}
	}
