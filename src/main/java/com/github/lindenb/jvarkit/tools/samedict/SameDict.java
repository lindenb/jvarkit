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
package com.github.lindenb.jvarkit.tools.samedict;

import java.io.BufferedReader;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.function.Function;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SequenceUtil.SequenceListsDifferException;

/**
BEGIN_DOC

## Return value

return -1 if there is no argument on the command line, if a dictionary cannot
be extracted from BAM/SAM/CRAM/BCF/VCF/FASTA etc...

return 0 if all argument share the same dictionary


## Example

```
$ find DIR -type f -name "*.bam" | java -jar dist/jvarkit.jar samedict


$ java -jar dist/jvarkit.jar samedict src/test/resources/S*.bam && echo "OK"
OK

$ java -jar dist/jvarkit.jar samedict || echo "ERROR"
[INFO][Launcher]samedict Exited with failure (-1)
ERROR

$ java -jar dist/jvarkit.jar samedict src/test/resources/S*.bam src/test/resources/ENCFF331CGL.rnaseq.b38.bam || echo "ERROR"
[INFO][Launcher]samedict Exited with failure (-1)
ERROR

$ java -jar dist/jvarkit.jar samedict src/test/resources/S*.bam src/test/resources/rotavirus_rf.unifiedgenotyper.vcf.gz src/test/resources/rotavirus_rf.fa && echo OK
OK

```



END_DOC
 */
@Program(name="samedict",
	description="check if all HTS files share the same dictionary",
	keywords={"dict","bed","sam","bam","vcf"},
	creationDate="20240724",
	modificationDate="20260218",
	jvarkit_amalgamion =  true
	)
public class SameDict extends Launcher {
	private static Logger LOG=Logger.of(SameDict.class);

	@Parameter(names={"--verbose"},description="Be verbose")
	private boolean be_verbose = false;

	
	@Override
	public int doWork(final List<String> args0) {
		
			
		try {
			final List<String> args;
 			if(args0.isEmpty() || (args0.size()==1 && args0.get(0).equals("-"))) {
 				args = new ArrayList<>();
				try(BufferedReader br  = IOUtils.openStreamForBufferedReader(stdin())) {
					String line;
					while((line=br.readLine())!=null)  {
						if(StringUtils.isBlank(line) || line.startsWith("#")) {
							args.add(line);
							}
						}
					}
 				}
 			else if(args0.size()==1 && args0.get(0).endsWith(".list")) {
 				args = new ArrayList<>();
				try(BufferedReader br  = IOUtils.openPathForBufferedReading(Paths.get(args0.get(0)))) {
					String line;
					while((line=br.readLine())!=null)  {
						if(StringUtils.isBlank(line) || line.startsWith("#")) {
							args.add(line);
							}
						}
					}
 				}
			else
				{
				args = args0;
				}
			
			SAMSequenceDictionary dict=null;
			for(final String filename: args) {
				final Optional<SAMSequenceDictionary> optDict = SequenceDictionaryUtils.extractDictionary(Paths.get(filename));
				if(!optDict.isPresent()) {
					if(be_verbose) LOG.warn("Cannot extract dict from "+filename);
					return -1;
					}
				if(dict==null) {
					dict = optDict.get();
					}
				else {
					try{
						SequenceUtil.assertSequenceDictionariesEqual(dict, optDict.get());
						}
					catch(final SequenceListsDifferException err) {
						if(be_verbose) {
							LOG.error(err);
							final Function<SAMSequenceRecord, String> toString=SSR->{
								if(SSR==null) return ".";
								return SSR.getSequenceName()+" len:"+SSR.getSequenceLength();
								};
							LOG.error("Dict in  "+filename+" is different from the others");
							final int nSeq = Math.max(optDict.get().size(), dict.size());
							for(int tid=0;tid < nSeq;++tid) {
								final SAMSequenceRecord recA  = tid<optDict.get().size()?optDict.get().getSequence(tid):null;
								final SAMSequenceRecord recB  = tid< dict.size() ? dict.getSequence(tid):null;
								if(recA==null || recB==null || !recA.equals(recB)) {
									LOG.error("$"+(tid+1)+"\t"+toString.apply(recA)+"\t"+toString.apply(recB));
									break;
									}
								}
							}	
						return -1;
						}
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
		new SameDict().instanceMainWithExit(args);
	}

}
