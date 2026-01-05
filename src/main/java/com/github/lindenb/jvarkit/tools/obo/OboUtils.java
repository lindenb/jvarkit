/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.obo;

import java.io.BufferedReader;
import java.io.File;
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
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.obo.OBOParser;
import com.github.lindenb.jvarkit.obo.OBOntology;

import htsjdk.samtools.util.StringUtil;

/**
BEGIN_DOC

## Example

```
$ java -jar dist/oboutils.jar --obo "https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo" -A "HP:0001626" | head
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



END_DOC
 *
 */
@Program(
		name="oboutils",
		description="OBO Ontology Utils.",
		keywords={"obo","ontology"},
		creationDate="20230105",
		modificationDate="20230105",
		jvarkit_amalgamion = true,
		menu="Utilities"
		)
public class OboUtils
	extends Launcher
	{
	private static Logger LOG=Logger.of(OboUtils.class);
	private enum Action{dump_table};
		
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-action","--action"},description= "What shoud I do ? default is dump as table.")
	private Action action = Action.dump_table;
	@Parameter(names={"-obo","--obo","--ontology"},description="Ontology in OBO format")
	private String hpoURI = "https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo";
	
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
	
		
	public OboUtils() {
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
		new OboUtils().instanceMainWithExit(args);
		}
	}
