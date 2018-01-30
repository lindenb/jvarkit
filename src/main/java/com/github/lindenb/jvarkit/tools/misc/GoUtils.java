/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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

*/package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.go.GoTree;
import com.github.lindenb.jvarkit.util.go.GoTree.RelType;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;

/**
BEGIN_DOC

## Example

children of  GO:0005216 'ion channel activity' 

```
$ java -jar dist/goutils.jar -R is_a  -A 'GO:0005216' 

#ACN	NAME	DEFINITION
GO:1905030	voltage-gated ion channel activity involved in regulation of postsynaptic membrane potential	Any voltage-gated ion channel activity that is involved in regulation of postsynaptic membrane potential.
GO:1905057	voltage-gated calcium channel activity involved in regulation of postsynaptic cytosolic calcium levels	Any voltage-gated calcium channel activity that is involved in regulation of postsynaptic cytosolic calcium ion concentration.
GO:1905054	calcium-induced calcium release activity involved in regulation of presynaptic cytosolic calcium ion concentration	Any calcium-induced calcium release activity that is involved in regulation of presynaptic cytosolic calcium ion concentration.
GO:1905058	calcium-induced calcium release activity involved in regulation of postsynaptic cytosolic calcium ion concentration	Any calcium-induced calcium release activity that is involved in regulation of postsynaptic cytosolic calcium ion concentration.
GO:0016286	small conductance calcium-activated potassium channel activity	Enables the transmembrane transfer of potassium by a channel with a unit conductance of 2 to 20 picoSiemens that opens in response to stimulus by internal calcium ions. Small conductance calcium-activated potassium channels are more sensitive to calcium than are large conductance calcium-activated potassium channels. Transport by a channel involves catalysis of facilitated diffusion of a solute (by an energy-independent process) involving passage through a transmembrane aqueous pore or channel, without evidence for a carrier-mediated mechanism.
GO:0043855	cyclic nucleotide-gated ion channel activity	Enables the transmembrane transfer of an ion by a channel that opens when a cyclic nucleotide has been bound by the channel complex or one of its constituent parts.
GO:0043854	cyclic nucleotide-gated mechanosensitive ion channel activity	Enables the transmembrane transfer of an ion by a channel that opens in response to a mechanical stress and when a cyclic nucleotide has been bound by the channel complex or one of its constituent parts.
GO:0099142	intracellularly ATP-gated ion channel activity	Enables the transmembrane transfer of an ion by a channel that opens when ATP has been bound by the channel complex or one of its constituent parts on the intracellular side of the plasma membrane.
GO:0099101	G-protein gated potassium channel activity	A potassium channel activity that is gated by binding of a G-protein beta-gamma dimer.
(...)
```

## Example

Use GO annotation to retrieve genes associated to GO:0005216 'ion channel activity' 

```
join -t $'\t' -1 1 -2 2 \
	<(java -jar dist/goutils.jar -A 'GO:0005216' | cut -f 1 | sort | uniq) \
	<(wget -q -O - "http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.goa_human.gz?rev=HEAD" | gunzip -c | grep -v '^!' | cut -f3,5 | uniq | LC_ALL=C sort -t $'\t' -k2,2) |\
sort -t $'\t' -k2,2 |\
grep SCN5A -A 10 -B 10
```

```
(...)
GO:0086006	SCN2B
GO:0005244	SCN3A
GO:0005248	SCN3A
GO:0005248	SCN3A
GO:0005248	SCN3B
GO:0086006	SCN3B
GO:0005248	SCN4A
GO:0005248	SCN4A
GO:0005248	SCN4B
GO:0086006	SCN4B
GO:0005244	SCN5A
GO:0005248	SCN5A
GO:0005248	SCN5A
GO:0005248	SCN5A
GO:0005248	SCN5A
GO:0005248	SCN5A
GO:0005248	SCN5A
GO:0086006	SCN5A
GO:0086060	SCN5A
GO:0086061	SCN5A
GO:0086062	SCN5A
GO:0086063	SCN5A
GO:0005248	SCN7A
GO:0005248	SCN7A
GO:0005248	SCN7A
GO:0005248	SCN7A
GO:0005248	SCN8A
GO:0005248	SCN8A
GO:0005248	SCN9A
GO:0005248	SCN9A
GO:0005248	SCN9A
GO:0005272	SCNN1A
(...)
```



END_DOC
 *
 */
@Program(
		name="goutils",
		description="Gene Ontology Utils. Retrieves terms from Gene Ontology",
		keywords={"geneontology","go"}
		)
public class GoUtils
	extends Launcher
	{
	private static Logger LOG=Logger.build(GoUtils.class).make(); 
		
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names= {"-G","--go","--gourl"},description=GoTree.GO_URL_OPT_DESC)
	private String goRdfUri = GoTree.GO_RDF_URL;
	@Parameter(names= {"-A","--accession",},description="User Go Terms accession numbers or name")
	private Set<String> userAccStrings = new HashSet<>();
	@Parameter(names= {"-R","--rel",},description="use those GO relationships. Default: all")
	private Set<String> relationTypesStr = new HashSet<>();
	@Parameter(names= {"-i","--inverse",},description="inverse the result")
	private boolean inverse=false;
    
	public GoUtils() {
	}
	
	@Override
	public int doWork(final List<String> args) {			
		PrintWriter out=null;
		try
			{
			if(!args.isEmpty())
				{
				LOG.error("too many arguments");
				return -1;
				}
			final Set<GoTree.RelType> relTypes;
			if(this.relationTypesStr.isEmpty())
				{
				relTypes=new HashSet<>(Arrays.asList(RelType.values()));
				}
			else
				{
				relTypes = this.relationTypesStr.stream().
						map(S->GoTree.RelType.valueOf(S)).
						collect(Collectors.toSet());
				}
			LOG.info("using rels:"+relTypes);
			LOG.info("parsing "+this.goRdfUri);
			final GoTree gotree=new GoTree.Parser().
					setDebug(false).
					setRelations(relTypes).
					parse(this.goRdfUri);
			final Set<GoTree.Term> userTerms = new HashSet<>();
			for(final String s:this.userAccStrings)
				{
				if(StringUtil.isBlank(s)) continue;
				GoTree.Term t= gotree.getTermByAccession(s);
				if(t==null)
					{
					t= gotree.getTermByName(s);
					if(t==null)
						{		
						LOG.error("cannot find user term "+s);
						return -1;
						}
					}
				userTerms.add(t);
				}
			out = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			out.println("#ACN\tNAME\tDEFINITION");
			for(final GoTree.Term t:gotree.getTerms())
				{
				boolean keep=false;
				
				if(userTerms.isEmpty())
					{
					keep=true;
					}
				else
					{
					for(final GoTree.Term userTerm:userTerms)
						{
						if(t.isDescendantOf(userTerm))
							{
							keep=true;
							break;
							}
						}
					}
				
				if(this.inverse) keep=!keep;
				if(keep)
					{
					out.print(t.getAcn());
					out.print('\t');
					out.print(t.getName());
					out.print('\t');
					out.print(t.getDefinition());
					out.println();
					}
				}
			
			out.flush();
			out.close();
			
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			}
		}
	
		
	public static void main(final String[] args)
		{
		new GoUtils().instanceMainWithExit(args);
		}
	}
