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
package com.github.lindenb.jvarkit.tools.ukbiobank;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/*
BEGIN_DOC

## Example

```
$ java -jar dist/ukbiobanksamples.jar  --ontology /path/to/coding19.tsv --tab  /path/to/ukb46178.tab --column "f.41270." -A 'I34'

[INFO][UKBiobankSelectSamples]using I349{I34.9 Nonrheumatic mitral valve disorder, unspecified} I34{I34 Nonrheumatic mitral valve disorders} I340{I34.0 Mitral (valve) insufficiency} I341{I34.1 Mitral (valve) prolapse} I342{I34.2 Nonrheumatic mitral (valve) stenosis} I348{I34.8 Other nonrheumatic mitral valve disorders}

1000767	R55 [R55 Syncope and collapse ]	W190 [W19.0 Home ]	Z961 [Z96.1 Presence of intraocular lens ]	K529 [K52.9 Non-infective gastro-enteritis and colitis, unspecified ]	K638 [K63.8 Other specified diseases of intestine ]	H041 [H04.1 Other disorders of lachrymal gland ]	I739 [I73.9 Peripheral vascular disease, unspecified ]	H269 [H26.9 Cataract, unspecified ]	R104 [R10.4 Other and unspecified abdominal pain ]	N946 [N94.6 Dysmenorrhoea, unspecified ]	N905 [N90.5 Atrophy of vulva ]	E119 [E11.9 Without complications ]	N908 [N90.8 Other specified noninflammatory disorders of vulva and perineum ]	N920 [N92.0 Excessive and frequent menstruation with regular cycle ]	I341 [I34.1 Mitral (valve) prolapse ] (...)
1002786	I10 [I10 Essential (primary) hypertension ]	R11 [R11 Nausea and vomiting ]	R12 [R12 Heartburn ]	K20 [K20 Oesophagitis ]	Z864 [Z86.4 Personal history of psychoactive substance abuse ]	M159 [M15.9 Polyarthrosis, unspecified ]	Z423 [Z42.3 Follow-up care involving plastic surgery of upper extremity ]	W010 [W01.0 Home ]	I120 [I12.0 Hypertensive renal disease with renal failure ]	J459 [J45.9 Asthma, unspecified ]	K299 [K29.9 Gastroduodenitis, unspecified ]	D509 [D50.9 Iron deficiency anaemia, unspecified ]	D125 [D12.5 Sigmoid colon ]	E835 [E83.5 Disorders of calcium metabolism ]	K298 [K29.8 Duodenitis ]	M170 [M17.0 Primary gonarthrosis, bilateral ]	I340 [I34.0 Mitral (valve) insufficiency ]...
1003404	R02 [R02 Gangrene, not elsewhere classified ]	E86 [E86 Volume depletion ]	J40 [J40 Bronchitis, not specified as acute or chronic ]	Z602 [Z60.2 Living alone ]I420 [I42.0 Dilated cardiomyopathy ]	I071 [I07.1 Tricuspid insufficiency ]	I340 [I34.0 Mitral (valve) insufficiency ]
1003703	M545 [M54.5 Low back pain ]	M754 [M75.4 Impingement syndrome of shoulder ]	I340 [I34.0 Mitral (valve) insufficiency ](...)
```

END_DOC
*/


@Program(name="ukbiobanksamples",
description="Select samples from ukbiobank",
keywords={"ukbiobank"},
creationDate="20210705",
modificationDate="20210706"
)
public class UKBiobankSelectSamples extends Launcher {
	private static final Logger LOG = Logger.build(UKBiobankSelectSamples.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--ontology","--coding"},description="coding/ontology file.",required=true)
	private Path ontologyPath = null;
	@Parameter(names={"--tab"},description="*.tab file.",required=true)
	private Path tabFilePath = null;
	@Parameter(names={"--column"},description="column pattern. Look for columns starts with 'x'",required=false)
	private String colPattern = "f.41270.";
	@Parameter(names={"--user-coding","-A"},description="limit to those coding. multiple separated by commas")
	private List<String> userCoding = new ArrayList<>();
	@Parameter(names={"--inverse"},description="inverse selection of sample when using option -A. get samples NOT having the phenotypes")
	private boolean inverse_selection = false;
	@Parameter(names={"--hide-phenotypes"},description="do not print phenotypes")
	private boolean hide_phenotypes = false;
	
	private String unquote(String s) {
		if(s.length()>1 && s.startsWith("\"") && s.endsWith("\"")) {
			return s.substring(1,s.length()-1);
		} else
		{
			return s;
		}
	}
	
	@Override
	public int doWork(List<String> args) {
		try {
			final UKBiobankOntology ontology = UKBiobankOntology.load(this.ontologyPath);
			final Set<UKBiobankOntology.Term> userTermsCoding = new HashSet<>();
			
			for(String userCoding : this.userCoding.stream().
				flatMap(S->Arrays.stream(CharSplitter.COMMA.split(S))).
				map(S->S.trim()).
				filter(S->!S.isEmpty()).
				collect(Collectors.toSet()))
				{
				final UKBiobankOntology.Term t = ontology.findTermByCoding(userCoding);
				if(t==null) {
					LOG.error("cannot fint user's coding "+userCoding+" in "+this.ontologyPath);
					return -1;
					}
				userTermsCoding.addAll(t.getAllChildrenIncludingSelf());
				}
			
			if(!userTermsCoding.isEmpty()) {
				LOG.info("using "+userTermsCoding.stream().map(S->S.toString()).collect(Collectors.joining(" ")));
			}
			
			try(BufferedReader br = IOUtils.openPathForBufferedReading(this.tabFilePath)) {
				String line = br.readLine();
				if(line==null) throw new IOException("cannot read first line of "+this.tabFilePath);
				String[] tokens = CharSplitter.TAB.split(line);
				if(!tokens[0].equals("f.eid")) throw new IOException("expected line to start with f.eid but got "+tokens[0]+" in "+this.tabFilePath);
				final List<Integer> columns_indexes = new ArrayList<>(tokens.length);
				for(int x=1;x< tokens.length;++x) {
					String label = tokens[x];
					if(!label.startsWith(this.colPattern)) continue;
					columns_indexes.add(x);
					}
				if(columns_indexes.isEmpty()) {
					LOG.error("found no columns starting with "+this.colPattern);
					return -1;
					}
				
				try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {

					while((line=br.readLine())!=null) {
						 tokens = CharSplitter.TAB.split(line);
						 //System.err.println(tokens[0]); 
						final Set<UKBiobankOntology.Term> terms = new HashSet<>();
						 for(Integer colidx: columns_indexes) {
							if(colidx >= tokens.length ) {
								throw new JvarkitException.TokenErrors(colidx.intValue()+1, tokens);
								}
							final String cell = unquote(tokens[colidx]);
							if(cell.equals("NA")) continue;
							final UKBiobankOntology.Term term = ontology.findTermByCoding(cell);
							if(term==null) {
								LOG.error("cannot find coding "+cell+" in "+this.ontologyPath+" for "+tokens[0]);
								return -1;
								}
							//LOG.info(term);
							terms.addAll(term.getAllChildrenIncludingSelf());
						 	//System.err.println("terms0 "+terms); 
							}
						//System.err.println(terms);
						boolean accept;
						
						if(userTermsCoding.isEmpty()) {
							accept =true;
							}
						else 
							{
							accept = terms.stream().anyMatch(T->userTermsCoding.contains(T));
							if(this.inverse_selection) accept = !accept;
							}
						//System.err.println(accept);
						if(!accept) continue;
						pw.print(tokens[0]);
						if(!hide_phenotypes) {
							//System.err.println("sort..");
							final List<UKBiobankOntology.Term> sortedterms = terms.stream().sorted((A,B)->Integer.compare(A.getDepth(), B.getDepth())).collect(Collectors.toList());
							//System.err.println("done sort");
							for(UKBiobankOntology.Term t:sortedterms) {
								pw.print("\t");
								pw.print(t.getCoding());
								pw.print(" [");
								pw.print(t.getMeaning());
								pw.print(" ]");
								}
							}
						pw.println();
						}//end while read
					pw.flush();
					}
				}
				
			
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new UKBiobankSelectSamples().instanceMainWithExit(args);

	}

}
