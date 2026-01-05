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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLine;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
/**
BEGIN_DOC

## EXAMPLE

```
$ cat input.bed | java -jar dist/addlinearindextobed.jar -R  human_g1k_v37.fasta 

10000   1       10000   177417
227417  1       227417  267719
317719  1       317719  471368
(...)
3060255274      Y       23951428        28819361
3095123207      Y       58819361        58917656
3095271502      Y       58967656        59363566
```

END_DOC
 */
@Program(
		name="addlinearindextobed",
		description="Use a Sequence dictionary to create a linear index for a BED file. Can be used as a X-Axis for a chart.",
		keywords={"bed","reference"},
		creationDate="20140201",
		modificationDate="20250115",
		jvarkit_amalgamion = true,
		menu="Deprecated/barely used"
		)
public class AddLinearIndexToBed extends Launcher
{
	private static final Logger LOG = Logger.of(AddLinearIndexToBed.class);

	@Parameter(names = { "-o", "--out" }, description = OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names = { "-R", "--reference","--dict"}, description = DICTIONARY_SOURCE, required = true)
	private Path refFile = null;
	@Parameter(names={"--regex"},description="keep chromosomes matching that regular expression. (use --ignore too prevent error the other chromosomes)")
	private String contig_regex="(chr)?[0-9XY]+";	
	@Parameter(names={"--ignore"},description="skip unknown chromosomes.")
	private boolean ignore_unknown_chromosomes=false;


	protected int doWork(BufferedReader in, PrintWriter out,final long tid2offset[],final SAMSequenceDictionary dictionary ) throws IOException {
		final CharSplitter tab = CharSplitter.TAB;
		final ContigNameConverter ctgNameConverter = ContigNameConverter.fromOneDictionary(dictionary);
		String line = null;
		while ((line = in.readLine()) != null) {
			if (BedLine.isBedHeader(line)) {
				out.println(line);
				continue;
			}
			if (line.isEmpty() || line.startsWith("#") || line.startsWith("track") || line.startsWith("browser"))
				continue;
			final String tokens[] = tab.split(line, 3);
			if (tokens.length < 2) {
				LOG.warn("Bad chrom/pos line:" + line);
				continue;
			}
			final String ctg =  ctgNameConverter.apply(tokens[0]);
			if(StringUtils.isBlank(ctg)) {
				throw new JvarkitException.ContigNotFoundInDictionary(tokens[0],dictionary);
				}
			
			final SAMSequenceRecord ssr = dictionary.getSequence(ctg);
			if (ssr == null) {
				if(this.ignore_unknown_chromosomes) continue;
				throw new JvarkitException.ContigNotFoundInDictionary(tokens[0],dictionary);
				}
			int pos0 = Integer.parseInt(tokens[1]);
			if (pos0 < 0 || pos0 >= ssr.getSequenceLength()) {
				LOG.warn("position is out of range for : " + line + " length(" + tokens[0] + ")="
						+ ssr.getSequenceLength());
			}
			out.print(tid2offset[ssr.getSequenceIndex()] + pos0);
			out.print('\t');
			out.print(line);
			out.println();
		}
		return 0;
	}

	@Override
	public int doWork(final List<String> args) {
		if (refFile == null) {
			LOG.error("Reference file undefined");
			return -1;
		}
		try {
			final SAMSequenceDictionary dict0 = SequenceDictionaryUtils.extractRequired(this.refFile);

			final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(
					dict0.getSequences().stream().
					filter(SSR->SSR.getContig().matches(contig_regex)).
					collect(Collectors.toList())
					);
			
			if(dictionary.isEmpty()) {
				LOG.error("empty dictionary");
				return -1;
				}
			
			final long[] tid2offset = new long[dictionary.size()];
			Arrays.fill(tid2offset, 0L);
			for (int i = 1; i < dictionary.size(); ++i) {
				tid2offset[i] = tid2offset[i - 1] + dictionary.getSequence(i - 1).getSequenceLength();
				}
			try(PrintWriter out = openPathOrStdoutAsPrintWriter(this.outputFile)) {
				if (args.isEmpty()) {
					try(BufferedReader br = IOUtils.openStreamForBufferedReader(stdin())) {
						doWork(br, out,tid2offset,dictionary);
					}
				} else {
					for (final String arg : args) {
						try(BufferedReader br = IOUtils.openURIForBufferedReading(arg)) {
							doWork(br, out,tid2offset,dictionary);
						}
					}
				}
				out.flush();
				}

			return 0;
		} catch (final Throwable err) {
			LOG.error(err);
			return -1;
		}
	}

	
	public static void main(final String[] args) {
		new AddLinearIndexToBed().instanceMainWithExit(args);
	}
}
