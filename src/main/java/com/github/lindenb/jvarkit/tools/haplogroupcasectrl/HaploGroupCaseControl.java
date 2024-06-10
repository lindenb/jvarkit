/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.haplogroupcasectrl;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**

BEGIN_DOC

### Examples

```
 java -jar jvarkit.jar haplogroupcasectrl --cases input.cases  --controls input.ctrls  --tree input.xml   input.tsv

```


END_DOC
*/
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.stats.FisherCasesControls;
import com.github.lindenb.jvarkit.pedigree.CasesControls;
import com.github.lindenb.jvarkit.phylotree.PhyloTree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
BEGIN_DOC

### Examples



END_DOC
 */
@Program(name="haplogroupcasectrl",
	description="Run Fisher test for Haplogroup input.",
	keywords={"haplogroup","burden","mitochondrial",""},
	modificationDate="20240610",
	creationDate="20240610",
	jvarkit_amalgamion = true
	)
public class HaploGroupCaseControl extends Launcher {
	private static final Logger LOG = Logger.build(HaploGroupCaseControl.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputPath = null;
	@Parameter(names={"--phylotree","--tree"},description=PhyloTree.OPT_DESC)
	private String phyloTreeUri = PhyloTree.DEFAULT_URI;
	@ParametersDelegate
	private CasesControls caseControls = new CasesControls();
	
	@Override
	public int doWork(List<String> args) {
		try {
			final PhyloTree phyloTree = PhyloTree.load(this.phyloTreeUri);
			this.caseControls.load();
			final Map<String,PhyloTree.HaploGroup> sample2haplogroup = new HashMap<>();
			try(BufferedReader br = super.openBufferedReader(oneFileOrNull(args))) {
				String line = br.readLine();
				if(line==null) throw new IOException("cannot read first line");
				final FileHeader header = new FileHeader(line, CharSplitter.TAB);
				header.assertColumnExists("sample");
				header.assertColumnExists("haplogroup");
				while((line=br.readLine())!=null) {
					final FileHeader.RowMap row = header.toMap(line);
					final String sample = row.get("sample");
					String hg = row.get("haplogroup");
					if(hg.contains("*")) hg=StringUtils.substringBefore(hg, "*");
					if(!phyloTree.hasHaplogroup(hg)) {
						LOG.error("haplogroup "+hg+" undefined in "+this.phyloTreeUri+". skipping");
						continue;
						}
					if(!this.caseControls.contains(sample)) {
						continue;
						}
					final PhyloTree.HaploGroup haplogroup  = phyloTree.getHaplogroupByName(hg);
					sample2haplogroup.put(sample, haplogroup);
					}
				}
			if(sample2haplogroup.isEmpty()) {
				LOG.error("no sample found / no intersection with cases/controls");
				return -1;
				}
			this.caseControls.retain(sample2haplogroup.keySet());
			this.caseControls.checkHaveCasesControls();
			
			try(PrintWriter w = super.openPathOrStdoutAsPrintWriter(outputPath)) {
				w.print("haplogroup");
				w.print("\t");
				w.print("n_cases_alt");
				w.print("\t");
				w.print("n_cases_ref");
				w.print("\t");
				w.print("n_ctrls_alt");
				w.print("\t");
				w.print("n_ctrls_ref");
				w.print("\t");
				w.print("fisher");
				w.println();
				
				for(PhyloTree.HaploGroup hg: phyloTree.getAllHaploGroups()) {
					final FisherCasesControls fisher = new FisherCasesControls(this.caseControls);
					for(final String sn: sample2haplogroup.keySet()) {
						final PhyloTree.HaploGroup hg_sn = sample2haplogroup.get(sn);
						if(hg.hasDescendant(hg_sn)) {
							fisher.accept(sn);
							}
						}
					
					w.print(hg.getName());
					w.print("\t");
					w.print(fisher.getCasesAltCount());
					w.print("\t");
					w.print(fisher.getCasesRefCount());
					w.print("\t");
					w.print(fisher.getControlsAltCount());
					w.print("\t");
					w.print(fisher.getControlsRefCount());
					w.print("\t");
					w.print(fisher.getAsDouble());
					w.print("\t");
					w.println();
					}
				
				w.flush();
				}
			return 0;
			} 
		catch (final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
		
	public static void main(String[] args) {
		new HaploGroupCaseControl().instanceMainWithExit(args);
		}
}
