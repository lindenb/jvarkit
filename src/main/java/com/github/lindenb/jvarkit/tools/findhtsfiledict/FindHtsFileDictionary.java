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
package com.github.lindenb.jvarkit.tools.findhtsfiledict;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;

/**
BEGIN_DOC

## Example


```
$ cat references.txt
name	fasta
rf	/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa
toy	/home/lindenb/src/jvarkit-git/src/test/resources/toy.fa


$ find . -name "*.bam" -o -name "*.vcf" -o -name "*.vcf.gz" -o -name "*.cram" |\
	java -jar dist/jvarkit.jar findhtsfiledict  -R references.txt

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
	java -jar dist/jvarkit.jar findhtsfiledict  -R references.txt -e ignore

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
	description="Scan a set of HTS files (VCF, BAM, CRAM, BCF, etc...), return a tab delimited file (path-of-file,path/url-to-fasta)",
	keywords={"sam","bam","cram","vcf","dict"},
	creationDate="20190912",
	modificationDate="20240824",
	jvarkit_amalgamion = true
	)
public class FindHtsFileDictionary extends Launcher {
	private static final Logger LOG = Logger.build(FindHtsFileDictionary.class).make();

	private enum OnError { ignore,fatal,empty};
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-D","-R","--references","--dictionaries"},description="A tab delimited file with a required header, and required columns: 'name' (name of the reference) and 'fasta' (path/url to the indexed fasta ref). An Optional column 'dict' for the path/url to the dict file",
			required = true)
	private Path repositoryFile = null;
	@Parameter(names={"-e","--errors"},description="What shall we do on error")
	private OnError onError = OnError.fatal;
	@Parameter(names={"--no-header"},description="Do not print header")
	private boolean disable_header=false;

	private final SequenceDictionaryExtractor dictExtractor=new SequenceDictionaryExtractor();

	
	private static class DictEntry
		{
		String name;
		SAMSequenceDictionary dict;
		String pathOrUrl;
		String dictPathOrUrl;
		}

	
 	
	private  List<DictEntry> readDictionaryRepository(final Path path) throws IOException {
		final List<DictEntry>  dictionaries = new ArrayList<>();
		IOUtil.assertFileIsReadable(path);
		FileHeader fileHeader;
		try(BufferedReader br=IOUtils.openPathForBufferedReading(path)) {
			String line= br.readLine();
			if(line==null) throw new IOException("cannot read first line of "+path);
			fileHeader=new FileHeader(line, S->CharSplitter.TAB.splitAsStringList(S));
			fileHeader.assertColumnExists("name");
			fileHeader.assertColumnExists("fasta");
			while((line=br.readLine())!=null) {
				final FileHeader.RowMap row = fileHeader.toMap(line);	
				final DictEntry entry = new DictEntry();
				entry.name = row.get("name");
				if(StringUtils.isBlank(entry.name)) {
					throw new IOException(" key \"name\" missing or empty for "+row);
					}
				if(dictionaries.stream().anyMatch(D->D.name.equalsIgnoreCase(entry.name))) {
					throw new IOException("duplicate dictionary named "+entry.name);
					}
						
				entry.pathOrUrl = row.get("fasta");
				if(dictionaries.stream().anyMatch(D->D.pathOrUrl.equals(entry.pathOrUrl))) {
					throw new IOException("duplicate fasta "+entry.pathOrUrl);
					}
				
						
				if(!IOUtil.isUrl(entry.pathOrUrl)) {
					final Path fasta=Paths.get(entry.pathOrUrl);
					IOUtil.assertFileIsReadable(fasta);
					
					if(!fasta.isAbsolute()) {
						LOG.warn("path is not absolute "+fasta);
						}
					}
				
				entry.dictPathOrUrl = row.getOrDefault("dict","");
				if(!StringUtils.isBlank(entry.dictPathOrUrl)) {
					if(dictionaries.stream().anyMatch(D->D.dictPathOrUrl.equals(entry.dictPathOrUrl))) {
						throw new IOException("duplicate dict "+entry.dictPathOrUrl);
						}
					
							
					if(!IOUtil.isUrl(entry.dictPathOrUrl)) {
						final Path dict =Paths.get(entry.dictPathOrUrl);
						IOUtil.assertFileIsReadable(dict);
						
						if(!dict.isAbsolute()) {
							LOG.warn("path is not absolute "+dict);
							}
						}
					entry.dict =dictExtractor.extractRequiredDictionary(entry.dictPathOrUrl);
					}
				else
					{
					// throw error if it's not supported fasta
					entry.dict =dictExtractor.extractRequiredDictionary(entry.pathOrUrl);
					}
				dictionaries.add(entry);
				}
			}
		return dictionaries;
		}
	
	
 	private void scan(final PrintWriter out,final String filename,final List<DictEntry> entries)
 		{ 		 		
 		final SAMSequenceDictionary dict;
 		try {
 			dict = this.dictExtractor.extractRequiredDictionary(filename);
 			} 
 		catch (final Throwable e)
			{
 			switch(onError)
 				{
 				case empty: out.println(String.join("\t", filename,".","."));break;
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
 			out.println(String.join("\t", filename ,entry.get().name,entry.get().pathOrUrl));
 			}
 		else
 			{
 			switch(onError)
				{
				case empty: out.println(String.join("\t", filename ,".","."));break;
				case ignore: break;
				default: throw new RuntimeIOException("No Reference found for "+filename);
				}
 			}
 		}
	
 	@Override
 	public int doWork(final List<String> args) {
		try
			{
			final List<DictEntry> entries = readDictionaryRepository(this.repositoryFile);
			try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)){
				if(!this.disable_header) out.println("#hts-file\tname\tfasta");
				if(args.isEmpty())
					{
					try(final BufferedReader in= IOUtils.openStdinForBufferedReader()) {
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
			return 0;
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
