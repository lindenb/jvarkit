/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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




import java.io.File;
import java.io.IOException;

import javax.script.Bindings;
import javax.script.CompiledScript;
import javax.script.ScriptException;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.fastq.FastqPairedWriter;
import com.github.lindenb.jvarkit.fastq.FastqRecordPair;
import com.github.lindenb.jvarkit.jcommander.OnePassFastqLauncher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.CloseableIterator;


/**
 * Pierre Lindenbaum PhD @yokofakun
 * @author lindenb
 * filters FASTQ with javascript


	BEGIN_DOC




### Motivation


The script puts 'rec' a FastqRecord, or 'pair' for an interleaved input, into the script context 



## Example

Find pairs of fastq where both reads contains a **PmeI** restriction site ( `GTTT/AAAC` )

```bash
$ paste <(gunzip -c F.fastq.gz | paste - - - -)  <(gunzip -c R.fastq.gz | paste - - - -) |\
  tr "\t" "\n" |\
  java -jar dist/fastqjs.jar -i -e 'pair.get(0).getReadString().contains("GTTTAAAC") && pair.get(1).getReadString().contains("GTTTAAAC") '

@HWI-1KL149:13:C0RNFACXX:8:1309:5373:60519 1:N:0:CTTGTA
TTCCAAAAATGTTTAAACTTTACAAATTTTCTTTCTGCAAAGGATATTTAAAACTTTGTCAAGACAAATATAAAAGTCTGTTCTTTTCATTAGTCTCTATA
+
CCCFFFFFHHHHHJJJJJJJJJBHIIJJJJJJJJJJIJIJJJJJJJJJJJJJJJJJJJJJJJIJJJJJJJJIJJIFEHIJHIJHHHHHHHFFFFFFFEEDE
@HWI-1KL149:13:C0RNFACXX:8:1309:5373:60519 2:N:0:CTTGTA
ACGCTTGATATTTGGTTTAAACATTTCTTGATTCAGAGAAGGTAGATGGTTATAGAGACTAATGAAAAGAACAGACTTTTATATTTGTCTTGACAAAGTTT
+
CCCFFFFFHHHHHJIIIJJIJJJJJIJJJJEHIJJJGIIIJJ?FHGIGIDDFHJJJJJJJIIJJIJJJJJJJJJIJHHHHHHHFFFFFFFEEEEEEDDDDC
@HWI-1KL149:13:C0RNFACXX:8:1309:6861:76085 1:N:0:CTTGTA
ACAGTATATCTATGTGAAAGTTAAAAAGAAATCGCTGTTTAGATGGAAGATGAGACCAGGTTATCATAGTTTTAGAAGAGGAGTTTAAACTTCATGCAGTG
+
CCCFDFFFHHHHHJIJJJJJIIJJJJJJJJJJJJJJJIJJJJIJJJJJJJJIIJJJJJJJFHIJJJJJJHIIJIJJHHHHHFFFFEEEEEEEDDEDDDDDD
@HWI-1KL149:13:C0RNFACXX:8:1309:6861:76085 2:N:0:CTTGTA
GGGTAGTCCACAAACAATGTGTTCATGTTGTCTCCCTCTTACTCACAAGCCTTCTGTAGCTCCCAACATTCACTGCATGAAGTTTAAACTCCTCTTCTAAA
+
@CCDDDFFHHHHHJJIJJJIJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJGIIJIJIIJJJJJJJJJIJJJIJJHHHHHHFFFFFFCEEEEEEDDDDDDED

```


END_DOC
	*/

@Program(name="fastqjs",
	keywords={"fastq"},
	description="Filters a FASTQ file using javascript( java nashorn engine). ",
	modificationDate="20220209"
	)
public class FastqJavascript
	extends OnePassFastqLauncher
	{
	private static final Logger LOG = Logger.build(FastqJavascript.class).make();

	@Parameter(names={"-N","--limit"},description="limit to 'N' records -1:all")
	private long LIMIT = -1L ;
	@Parameter(names={"-e"},description="javascript expression")
	private String javascriptExpr = null;
	@Parameter(names={"-f"},description="javascript file")
	private File javascriptFile = null;

	
	
	private CompiledScript  script=null;
	//private File failingReadsFile=null;
	//
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	/** mutable fastq record */
	public static class Record
		{
		private long nLine;
		private String name;
		private String sequence;
		private String name2;
		private String qualities;
		
		Record(final FastqRecord rec)
			{
			name=rec.getReadName();
			sequence=rec.getReadString();
			name2=rec.getBaseQualityHeader();
			qualities=rec.getBaseQualityString();
			}
		public String getReadHeader() { return name;}
		public String getReadString() { return sequence;}
		public String getBaseQualityHeader() { return name2;}
		public String getBaseQualityString() { return qualities;}
		public long getLine() {
			return nLine;
			}
		public void setName(String name) {
			this.name = name;
			}
		public void setReadString(String sequence) {
			this.sequence = sequence;
		}
		public void setBaseQualityHeader(String name2) {
			this.name2 = name2;
		}
		public void setBaseQualityString(String qualities) {
			this.qualities = qualities;
			}
		private FastqRecord toFastqRecord()
			{
			return new FastqRecord(name,sequence,name2,qualities);
			}
		@Override
		public int hashCode() {
			return toFastqRecord().hashCode();
			}
		@Override
		public String toString() {
			return toFastqRecord().toString();
			}
		}
	
	/** a pair of FASTQ */
	public static class Pair
		{
		private long nLine;
		private Record rec1;
		private Record rec2;
		Pair(FastqRecordPair pair)
			{
			this.rec1=new Record(pair.get(0));
			this.rec2=new Record(pair.get(1));
			}
		public long getLine() {
			return nLine;
			}
		public Record getFirst() { return get(0);}
		public Record getSecond() { return get(1);}
		public Record get(int index)
			{
			switch(index)
				{
				case 0: return rec1;
				case 1: return rec2;
				default: throw new IllegalArgumentException("index:="+index);
				}
			}
		
		@Override
		public int hashCode() {
			return get(0).hashCode()*31+get(1).hashCode();
			}
		
		@Override
		public String toString() {
			return get(0).toString()+"\n"+get(1).toString();
			}
		}
	
	

	
	@Override
	protected int beforeFastq() {
		try {
			this.script  = super.compileJavascript(this.javascriptExpr,this.javascriptFile);
			}
		catch(final Exception err) {
			getLogger().error(err);
			return -1;
			}
		

		return 0;
		}
	
	
	@Override
	protected int runSingleEnd(FastqReader iter, FastqWriter fws) throws IOException {
		try {
			long count=0L;
			final Bindings bindings = this.script.getEngine().createBindings();
			while(iter.hasNext()) {
				final Record record=new Record(iter.next());
				record.nLine=count;
				bindings.put("rec", record);
				if(super.evalJavaScriptBoolean(this.script, bindings))
					{
					fws.write(record.toFastqRecord());
					}
				++count;
				if(this.LIMIT>0L && count>=this.LIMIT) break;
				}
			return 0;
			}
		catch(final ScriptException err) {
			getLogger().error(err);
			return -1;
			}
		}
	
	
	@Override
	protected int runPairedEnd(CloseableIterator<FastqRecordPair> iter, FastqPairedWriter fws) throws IOException {
		try {
			long count=0L;
			final Bindings bindings = this.script.getEngine().createBindings();
			while(iter.hasNext()) {
				final Pair p = new Pair(iter.next());
				bindings.put("pair", p);
				if(super.evalJavaScriptBoolean(this.script, bindings))
					{
					fws.write(p.rec1.toFastqRecord(),p.rec2.toFastqRecord());
					}
				++count;
				if(this.LIMIT>0L && count>=this.LIMIT) break;
				}
			return 0;
			}
		catch(final ScriptException err) {
			getLogger().error(err);
			return -1;
			}
		}
			
			
		
	public static void main(String[] args) throws Exception
		{
		new FastqJavascript().instanceMainWithExit(args);
		}
	
	}
