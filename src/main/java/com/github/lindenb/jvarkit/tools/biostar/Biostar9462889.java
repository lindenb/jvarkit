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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.biostar;

import java.util.Arrays;
import java.util.Collections;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.AbstractBamSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMRecord;


/**

BEGIN_DOC

##Example

```
$ java -jar dist/biostar9462889.jar --manifest jeter.mf -o TMP --regex '^(RF[0-9]+)_' src/test/resources/S1.bam
WARNING: BAM index file /home/lindenb/src/jvarkit/src/test/resources/S1.bam.bai is older than BAM /home/lindenb/src/jvarkit/src/test/resources/S1.bam
[INFO][Biostar9462889]Creating output for "RF01" N=1
[INFO][Biostar9462889]Creating output for "RF02" N=2
[INFO][Biostar9462889]Creating output for "RF03" N=3
[INFO][Biostar9462889]Creating output for "RF04" N=4
[INFO][Biostar9462889]Creating output for "RF05" N=5
[INFO][Biostar9462889]Creating output for "RF06" N=6
[INFO][Biostar9462889]Creating output for "RF07" N=7
[INFO][Biostar9462889]Creating output for "RF08" N=8
[INFO][Biostar9462889]Creating output for "RF09" N=9
[INFO][Biostar9462889]Creating output for "RF10" N=10
[INFO][Biostar9462889]Creating output for "RF11" N=11
[WARN][Biostar9462889]0 read(s) where lost because the regex '^(RF[0-9]+)_' failed.

$ ls TMP/*.bam
TMP/split.000001.bam  TMP/split.000004.bam  TMP/split.000007.bam  TMP/split.000010.bam
TMP/split.000002.bam  TMP/split.000005.bam  TMP/split.000008.bam  TMP/split.000011.bam
TMP/split.000003.bam  TMP/split.000006.bam  TMP/split.000009.bam

$ cat jeter.mf
RF01	TMP/split.000001.bam
RF02	TMP/split.000002.bam
RF03	TMP/split.000003.bam
RF04	TMP/split.000004.bam
RF05	TMP/split.000005.bam
RF06	TMP/split.000006.bam
RF07	TMP/split.000007.bam
RF08	TMP/split.000008.bam
RF09	TMP/split.000009.bam
RF10	TMP/split.000010.bam
RF11	TMP/split.000011.bam
```

END_DOC
*/
@Program(name="biostar9462889",
	keywords={"sam","bam","split","util"},
	description="Extracting reads from a regular expression in a bam file",
	creationDate="20210402",
	modificationDate="20210402",
	biostars=9462889
	)
public class Biostar9462889 extends AbstractBamSplitter<Biostar9462889.Tokens>
	{
	private static final Logger LOG = Logger.build(Biostar9462889.class).make();

	@Parameter(names= {"--regex","-regex"},description="Regular expression that can be used to parse read names in the incoming SAM file. Regex groups are used to classify the reads.",required=true)
	private String regexStr="";
	
	private Pattern compiledRegex = null;
	
	static class Tokens {
		final String[] tokens;
		Tokens(final Matcher matcher) {
			this.tokens = new String[matcher.groupCount()];
			for(int i=0;i< matcher.groupCount();i++) {
				this.tokens[i] = matcher.group(i+1);/** yes i+1 that's weird but group(0) is whole string */
				}
			}
		@Override
		public int hashCode() {
			return Arrays.hashCode(tokens);
			}
		@Override
		public boolean equals(Object obj) {
			if(obj==this)return true;
			return Arrays.equals(this.tokens,Tokens.class.cast(obj).tokens);
			}
		@Override
		public String toString() {
			return String.join(" ", this.tokens);
			}
		}
	
	protected Logger getLogger() {
		return LOG;
	};
	
	@Override
	protected int beforeSam() {
		try {
			this.compiledRegex = Pattern.compile(this.regexStr);
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		return super.beforeSam();
		}
	
	@Override
	protected Set<Tokens> createKeys(final SAMRecord rec) {
		final String rn = rec.getReadName();
		final Matcher matcher = this.compiledRegex.matcher(rn);
		if(!matcher.find() || matcher.groupCount()==0) return Collections.emptySet();
		return Collections.singleton(new Tokens(matcher));
		}
	
	public static void main(final String[] args) {
		new Biostar9462889().instanceMainWithExit(args);
	}

}
