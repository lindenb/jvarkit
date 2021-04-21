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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

/**
BEGIN_DOC

## Example


```
$ cat references.txt

name : rf
fasta : /home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa

name : toy
fasta : /home/lindenb/src/jvarkit-git/src/test/resources/toy.fa


$ find . -name "*.bam" -o -name "*.vcf" -o -name "*.vcf.gz" -o -name "*.cram" |\
	java -jar dist/findhtsfiledict.jar  -R references.txt

[SEVERE][FindHtsFileDictionary]No Reference found for ./e90235dca2ff2d151ef796ad0be98915/test01.vcf
java.io.IOException: No Reference found for ./e90235dca2ff2d151ef796ad0be98915/test01.vcf
	at com.github.lindenb.jvarkit.tools.misc.FindHtsFileDictionary.scan(FindHtsFileDictionary.java:193)
	at com.github.lindenb.jvarkit.tools.misc.FindHtsFileDictionary.doWork(FindHtsFileDictionary.java:213)
	at com.github.lindenb.jvarkit.util.jcommander.Launcher.instanceMain(Launcher.java:756)
	at com.github.lindenb.jvarkit.util.jcommander.Launcher.instanceMainWithExit(Launcher.java:919)
	at com.github.lindenb.jvarkit.tools.misc.FindHtsFileDictionary.main(FindHtsFileDictionary.java:241)
#hts-file	dict.name	fasta
/home/lindenb/src/jvarkit-git/./b7/83f96c410c7cd75bc732d44a1522a7/Gene_32_1507.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
[INFO][Launcher]findhtsfiledict Exited with failure (-1)


$ find ${PWD}/ -name "*.bam" -o -name "*.vcf" -o -name "*.vcf.gz" -o -name "*.cram" |\
	java -jar dist/findhtsfiledict.jar  -R references.txt -e ignore

#hts-file	dict.name	fasta
/home/lindenb/src/jvarkit-git/b7/83f96c410c7cd75bc732d44a1522a7/Gene_32_1507.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/db/aee9cc8f5c9c3d39c7af4cec63b7a5/Gene_0_1061.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/jeter.cram	toy	/home/lindenb/src/jvarkit-git/src/test/resources/toy.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S5.bam	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/toy.vcf.gz	toy	/home/lindenb/src/jvarkit-git/src/test/resources/toy.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S1.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S2.bam	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S5.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S4.bam	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S3.bam	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/toy.cram	toy	/home/lindenb/src/jvarkit-git/src/test/resources/toy.fa
/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.unifiedgenotyper.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S2.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S4.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S1.bam	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/toy.bam	toy	/home/lindenb/src/jvarkit-git/src/test/resources/toy.fa
/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.ann.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/S3.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.freebayes.vcf.gz	rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
```



END_DOC
 */


@Program(name="findhtsfiledict",
	description="Scan a set of HTS files (VCF, BAM, CRAM), return a tab delimited file (path-of-file,path-to-fasta)",
	keywords={"sam","bam","cram","vcf","dict"},
	creationDate="20190912",
	modificationDate="20190912"
	)
public class FindHtsFileDictionary extends Launcher {
	private static final Logger LOG = Logger.build(FindHtsFileDictionary.class).make();

	private enum OnError { ignore,raise,empty};
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-D","-R","--repository","--dictionaries"},description="A repository is a set of record separated by a white space. "
			+ "Each record is a line of 'key : value' or 'key=value'. Two keys are required: name and fasta ",
			required = true)
	private Path repositoryFile = null;
	@Parameter(names={"-e","--errors"},description="What shall we do on error")
	private OnError onError = OnError.raise;

	private static class DictEntry
		{
		String name;
		SAMSequenceDictionary dict;
		Path fasta;
		}

	
 	
	private  List<DictEntry> readDictionaryRepository(final Path path) throws IOException {
		final List<DictEntry>  dictionaries = new ArrayList<>();
		final Map<String,String> properties= new HashMap<>();
		IOUtil.assertFileIsReadable(path);
		try(BufferedReader br=IOUtils.openPathForBufferedReading(path)) {
			String line;
			for(;;) {
				line= br.readLine();
				if(StringUtils.isBlank(line)) {
					if(!properties.isEmpty()) {
						final DictEntry entry = new DictEntry();
						entry.name = properties.get("name");
						if(StringUtils.isBlank(entry.name)) {
							throw new IOException(" key \"name\" missing or empry for "+properties);
							}
						if(dictionaries.stream().anyMatch(D->D.name.equalsIgnoreCase(entry.name))) {
							throw new IOException("duplicate dictionary named "+entry.name);
							}
						
						if(!properties.containsKey("fasta")) {
							throw new IOException(" key \"fasta\" missing for "+properties);
							}
						entry.fasta = Paths.get(properties.get("fasta"));
						if(dictionaries.stream().anyMatch(D->D.fasta.equals(entry.fasta))) {
							throw new IOException("duplicate fasta "+entry.fasta);
							}
						
						
						IOUtil.assertFileIsReadable(entry.fasta);
						
						if(!entry.fasta.isAbsolute()) {
							LOG.warn("path is not absolute "+entry.fasta);
							}
						// throw error if it's not supported fasta
						ReferenceSequenceFileFactory.getFastaExtension(entry.fasta);
						entry.dict =SequenceDictionaryUtils.extractRequired(entry.fasta);
						dictionaries.add(entry);
						}
					if(line==null) break;
					properties.clear();
					continue;
					}
				int sep = line.indexOf(':');
				if(sep==-1) sep=line.indexOf('=');
				if(sep==-1) throw new IOException("cannot find separator ':' or '=' in "+line);
				final String left = line.substring(0,sep).trim().toLowerCase();
				final String right = line.substring(sep+1).trim();
				if(properties.containsKey(left)) {
					throw new IOException("duplicate key "+left+" for "+properties);
					}
				properties.put(left, right);
				}
			}
		return dictionaries;
		}
	
	
 	private void scan(final PrintWriter out,final String filename,final List<DictEntry> entries)
 		{ 		
 		final Path f=Paths.get(filename);
 		
 		final SAMSequenceDictionary dict;
 		try {
 			dict = SAMSequenceDictionaryExtractor.extractDictionary(f);
 			} 
 		catch (final Throwable e)
			{
 			switch(onError)
 				{
 				case empty: out.println(String.join("\t", f.toAbsolutePath().toString(),".","."));break;
 				case ignore: break;
 				default: throw new RuntimeIOException(e);
 				}
 			return;
			}
 		final Optional<DictEntry> entry = entries.
 				stream().
 				filter(E->SequenceUtil.areSequenceDictionariesEqual(E.dict, dict)).
 				findFirst(); 
 		
 		if(entry.isPresent()) {
 			out.println(String.join("\t", f.toAbsolutePath().toString(),entry.get().name,entry.get().fasta.toAbsolutePath().toString()));
 			}
 		else
 			{
 			switch(onError)
				{
				case empty: out.println(String.join("\t", f.toAbsolutePath().toString(),".","."));break;
				case ignore: break;
				default: throw new RuntimeIOException("No Reference found for "+f);
				}
 			}
 		}
	
 	@Override
 	public int doWork(final List<String> args) {
		try
			{
			final List<DictEntry> entries = readDictionaryRepository(this.repositoryFile);
			try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)){
				out.println("#hts-file\tdict.name\tfasta");
				if(args.isEmpty())
					{
					try(final BufferedReader in=new BufferedReader(new InputStreamReader(stdin()))) {
						in.lines().
							filter(line->!(StringUtils.isBlank(line) || line.startsWith("#"))).
							forEach(line->scan(out,line,entries));						
						}
					}
				else
					{
					for(final String filename:args)
						{
						scan(out,filename,entries);
						}
					}
				out.flush();
				}
			return RETURN_OK;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		}

 	public static void main(String[] args)
		{
		new FindHtsFileDictionary().instanceMainWithExit(args);
		}

}
