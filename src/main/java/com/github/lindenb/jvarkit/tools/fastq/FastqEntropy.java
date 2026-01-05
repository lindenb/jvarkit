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
package com.github.lindenb.jvarkit.tools.fastq;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.zip.Deflater;
import java.util.zip.GZIPOutputStream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.fastq.FastqRecordPair;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.jcommander.FastqLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.Counter;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;

/**
BEGIN_DOC

## Example

```bash
$  java -jar dist/fastqentropy.jar file.fastq.gz
```
output (length/count):
```
36	1
37	1
38	2
39	2
128	3
41	7
40	9
127	10
42	18
43	20
44	34
126	35
45	71
46	114
125	114
47	198
124	254
48	272
49	475
123	552
50	579
51	787
52	1028
122	1085
53	1281
54	1391
55	1624
121	1846
56	1847
57	1907
58	2182
60	2234
59	2245
62	2387
61	2392
63	2460
64	2478
66	2613
65	2655
67	2740
68	2753
78	2881
77	2919
79	2921
76	2928
75	2954
80	2958
69	2969
70	3014
81	3014
71	3030
72	3056
82	3067
74	3073
73	3075
120	3157
83	3175
84	3381
85	3433
86	3627
87	3846
88	4098
89	4417
90	4535
119	4951
91	5096
92	5544
93	6045
94	6475
95	6913
118	7089
96	7211
97	7697
98	7877
99	8118
100	8194
101	8235
102	8280
103	8496
104	8740
105	9144
117	9335
106	9571
107	10767
116	11688
108	12174
109	13664
115	14461
110	15637
114	16983
111	17347
113	18332
112	18698
```
END_DOC

 */

@Program(name="fastqentropy",
	description="Compute the Entropy of a Fastq file (distribution of the length(gzipped(sequence))",
	keywords={"fastq"},
	modificationDate="20220209"
	)
public class FastqEntropy extends FastqLauncher
	{
	private static final Logger LOG = Logger.of(FastqEntropy.class);

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path fileout = null;

	private final Counter<Long> length2count=new Counter<Long>();
	
	private static class BestCompressionOutputStream extends GZIPOutputStream
		{
		BestCompressionOutputStream() throws IOException
			{
			super(new NullOuputStream());
			def.setLevel(Deflater.BEST_COMPRESSION);
			}
		public long getByteWrittenCount()
			{
			return NullOuputStream.class.cast(super.out).getByteWrittenCount();
			}
		}

	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	
	private void convert(FastqRecord rec)
		{
		try {
			try(BestCompressionOutputStream gzout=new BestCompressionOutputStream()) {
				gzout.write(rec.getBaseQualityString().getBytes());
				gzout.flush();
				this.length2count.incr(gzout.getByteWrittenCount());
				}
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	
	
	@Override
	protected int runPairedEnd(CloseableIterator<FastqRecordPair> iter) throws IOException {
		while(iter.hasNext()) {
			iter.next().forEach(REC->convert(REC));
			}
		return dump();
		}
	@Override
	protected int runSingleEnd(FastqReader iter) throws IOException {
		while(iter.hasNext()) {
			convert(iter.next());
			}
		return dump();
		}
	
	private int dump() throws IOException {
		try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.fileout)) {
			for(Long n:this.length2count.keySetIncreasing())
				{
				pw.print(n);
				pw.print('\t');
				pw.println(this.length2count.count(n));
				if(pw.checkError()) break;
				}
			pw.flush();
			}
		return 0;
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new FastqEntropy().instanceMainWithExit(args);
	}

}
