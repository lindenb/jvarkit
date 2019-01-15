/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
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
		keywords={"bed","reference"}
		)
public class AddLinearIndexToBed extends Launcher
{
	private static final Logger LOG = Logger.build(AddLinearIndexToBed.class).make();

	@Parameter(names = { "-o", "--out" }, description = OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names = { "-R", "--reference" }, description = INDEXED_FASTA_REFERENCE_DESCRIPTION, required = true)
	private File refFile = null;

	public AddLinearIndexToBed() {

	}

	private SAMSequenceDictionary dictionary = null;
	private long tid2offset[] = null;

	protected int doWork(InputStream is, PrintStream out) throws IOException {
		final Pattern tab = Pattern.compile("[\t]");
		BufferedReader in = new BufferedReader(new InputStreamReader(is));
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
			final SAMSequenceRecord ssr = this.dictionary.getSequence(tokens[0]);
			if (ssr == null) {
				for (SAMSequenceRecord sr2 : this.dictionary.getSequences()) {
					LOG.info("available " + sr2.getSequenceName());
				}
				throw new IOException("undefined chromosome:" + tokens[0]);
			}
			int pos0 = Integer.parseInt(tokens[1]);
			if (pos0 < 0 || pos0 >= ssr.getSequenceLength()) {
				LOG.warn("position is out of range for : " + line + " length(" + tokens[0] + ")="
						+ ssr.getSequenceLength());
			}
			out.print(this.tid2offset[ssr.getSequenceIndex()] + pos0);
			out.print('\t');
			out.print(line);
			out.println();
			if (out.checkError())
				break;
		}
		return 0;
	}

	@Override
	public int doWork(final List<String> args) {
		if (refFile == null) {
			LOG.error("Reference file undefined");
			return -1;
		}
		PrintStream out = null;
		try {
			htsjdk.samtools.reference.IndexedFastaSequenceFile ref = new IndexedFastaSequenceFile(this.refFile);
			this.dictionary = ref.getSequenceDictionary();
			ref.close();

			if (this.dictionary == null) {
				throw new JvarkitException.FastaDictionaryMissing(this.refFile);
			}

			this.tid2offset = new long[this.dictionary.size()];
			Arrays.fill(this.tid2offset, 0L);
			for (int i = 1; i < this.dictionary.size(); ++i) {
				this.tid2offset[i] = this.tid2offset[i - 1] + this.dictionary.getSequence(i - 1).getSequenceLength();
			}
			out = openFileOrStdoutAsPrintStream(this.outputFile);

			if (args.isEmpty()) {
				LOG.info("reading stdin");
				doWork(stdin(), out);
			} else {
				for (final String arg : args) {
					LOG.info("opening " + arg);
					InputStream in = IOUtils.openURIForReading(arg);
					doWork(in, out);
					CloserUtil.close(out);
				}
			}
			out.flush();

			return RETURN_OK;
		} catch (Exception err) {
			LOG.error(err);
			return -1;
		} finally {
			CloserUtil.close(out);
			dictionary = null;
			tid2offset = null;
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new AddLinearIndexToBed().instanceMainWithExit(args);
	}
}
