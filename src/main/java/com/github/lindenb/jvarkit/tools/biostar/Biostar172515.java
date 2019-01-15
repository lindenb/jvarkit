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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.OutputStream;
import java.util.List;
import java.util.Set;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.BAMFileSpan;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.Bin;
import htsjdk.samtools.BinList;
import htsjdk.samtools.BrowseableBAMIndex;
import htsjdk.samtools.Chunk;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;

/**

BEGIN_DOC

### Example

``` 
$ find DIR -name "*.bam" | xargs java -jar dist/biostar172515.jar  | xmllint --format -

<?xml version="1.0" encoding="UTF-8"?>
<bai-list>
<bam bam="DIR/exampleBAM.bam" has-index="true" n_ref="1">
    <reference ref-id="0" ref-name="chr1" ref-length="100000" n_aligned="33" n_bin="12" n_no_coor="0">
      <bin first-locus="1" last-locus="536870912" level="0" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="67108864" level="1" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="8388608" level="2" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="1048576" level="3" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="131072" level="4" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="16384" level="5" first-offset="828" n_chunk="1">
        <chunk chunk_beg="828" chunk_end="1963"/>
      </bin>
      <bin first-locus="16385" last-locus="32768" level="5" first-offset="1963" n_chunk="1">
        <chunk chunk_beg="1963" chunk_end="3323"/>
      </bin>
      <bin first-locus="32769" last-locus="49152" level="5" first-offset="3323" n_chunk="1">
        <chunk chunk_beg="3323" chunk_end="4687"/>
      </bin>
      <bin first-locus="49153" last-locus="65536" level="5" first-offset="4687" n_chunk="1">
        <chunk chunk_beg="4687" chunk_end="6501"/>
      </bin>
      <bin first-locus="65537" last-locus="81920" level="5" first-offset="0" n_chunk="0"/>
      <bin first-locus="81921" last-locus="98304" level="5" first-offset="6501" n_chunk="1">
        <chunk chunk_beg="6501" chunk_end="238223360"/>
      </bin>
      <bin first-locus="98305" last-locus="114688" level="5" first-offset="0" n_chunk="0"/>
    </reference>
  </bam>
  <bam bam="/DIR/exampleBAM2.bam" has-index="true" n_ref="1">
    <reference ref-id="0" ref-name="chr1" ref-length="100000" n_aligned="33" n_bin="12" n_no_coor="0">
      <bin first-locus="1" last-locus="536870912" level="0" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="67108864" level="1" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="8388608" level="2" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="1048576" level="3" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="131072" level="4" first-offset="0" n_chunk="0"/>
      <bin first-locus="1" last-locus="16384" level="5" first-offset="828" n_chunk="1">
        <chunk chunk_beg="828" chunk_end="1963"/>
      </bin>
      <bin first-locus="16385" last-locus="32768" level="5" first-offset="1963" n_chunk="1">
        <chunk chunk_beg="1963" chunk_end="3323"/>
      </bin>
      <bin first-locus="32769" last-locus="49152" level="5" first-offset="3323" n_chunk="1">
        <chunk chunk_beg="3323" chunk_end="4687"/>
      </bin>
      <bin first-locus="49153" last-locus="65536" level="5" first-offset="4687" n_chunk="1">
        <chunk chunk_beg="4687" chunk_end="6501"/>
      </bin>
      <bin first-locus="65537" last-locus="81920" level="5" first-offset="0" n_chunk="0"/>
      <bin first-locus="81921" last-locus="98304" level="5" first-offset="6501" n_chunk="1">
        <chunk chunk_beg="6501" chunk_end="238223360"/>
      </bin>
      <bin first-locus="98305" last-locus="114688" level="5" first-offset="0" n_chunk="0"/>
    </reference>
  </bam>
</bai-list>
```

END_DOC
*/


@Program(name="biostar172515",
description="Convert BAI to XML",
biostars=172515,
keywords={"bai","bam","xml"}
)
public class Biostar172515 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar172515.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	private XMLStreamWriter w=null;
	public Biostar172515() {
		}
	
	@Override
	public int doWork(final List<String> inputFiles) {
		final SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().
					setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, Boolean.TRUE).
					validationStringency(ValidationStringency.LENIENT);
		OutputStream stream=null;
		SamReader samReader = null;
		Set<String> args=IOUtils.unrollFiles(inputFiles);
		try {
			stream = super.openFileOrStdoutAsStream(this.outputFile);
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			this.w = xof.createXMLStreamWriter(stream);
			this.w.writeStartDocument("UTF-8", "1.0");
			this.w.writeStartElement("bai-list");
			for(final String filename:args)
				{
				this.w.writeStartElement("bam");
				this.w.writeAttribute("bam", filename);
				samReader = samReaderFactory.open(SamInputResource.of(filename));
				this.w.writeAttribute("has-index",String.valueOf(samReader.hasIndex()));
				
				if(!samReader.hasIndex())
					{
					this.w.writeEndElement();
					samReader.close();
					continue;
					}
				final SamReader.Indexing indexing = samReader.indexing();

				
				if(!indexing.hasBrowseableIndex())
					{
					this.w.writeComment("no browseable index");
					this.w.writeEndElement();
					samReader.close();
					continue;
					}
				
				final SAMSequenceDictionary dict= samReader.getFileHeader().getSequenceDictionary();
				this.w.writeAttribute("n_ref",String.valueOf(dict.size()));

				
				
				final BrowseableBAMIndex baiFile;
				try {
					baiFile = indexing.getBrowseableIndex();
					}
				catch(Exception err)
					{
					this.w.writeComment("no browseable index");
					this.w.writeEndElement();
					samReader.close();
					continue;
					}
				for(int tid=0;tid< dict.size();++tid)
					{
					final SAMSequenceRecord ssr = dict.getSequence(tid);
					final BAMIndexMetaData baiMetaData = baiFile.getMetaData(tid);

					
					this.w.writeStartElement("reference");
					this.w.writeAttribute("ref-id",String.valueOf(tid));
					this.w.writeAttribute("ref-name",ssr.getSequenceName());
					this.w.writeAttribute("ref-length",String.valueOf(ssr.getSequenceLength()));
					this.w.writeAttribute("n_aligned",String.valueOf(baiMetaData.getAlignedRecordCount()));
					BinList binList =baiFile.getBinsOverlapping(tid, 1, ssr.getSequenceLength());
					int n_bin=0;
					for(@SuppressWarnings("unused") final Bin binItem:binList) n_bin++;
					this.w.writeAttribute("n_bin",String.valueOf(n_bin));
					
					
					this.w.writeAttribute("n_no_coor",String.valueOf(baiMetaData.getUnalignedRecordCount()));
					for(final Bin binItem:binList)
						{
						
						this.w.writeStartElement("bin");
						this.w.writeAttribute("first-locus",String.valueOf(baiFile.getFirstLocusInBin(binItem)));
						this.w.writeAttribute("last-locus",String.valueOf(baiFile.getLastLocusInBin(binItem)));
						this.w.writeAttribute("level",String.valueOf(baiFile.getLevelForBin(binItem)));
						final BAMFileSpan span= baiFile.getSpanOverlapping(binItem);
						this.w.writeAttribute("first-offset",String.valueOf(span.getFirstOffset()));
						final List<Chunk> chunks = span.getChunks();
						this.w.writeAttribute("n_chunk",String.valueOf(chunks.size()));

						
						for(final Chunk chunk:chunks)
							{
							this.w.writeEmptyElement("chunk");
							this.w.writeAttribute("chunk_beg",String.valueOf(chunk.getChunkStart()));
							this.w.writeAttribute("chunk_end",String.valueOf(chunk.getChunkEnd()));
							}
						
						this.w.writeEndElement();
						}
										
					this.w.writeEndElement();
					}
				
				this.w.writeEndElement();
				samReader.close();
				}
			this.w.writeEndElement();
			this.w.flush();
			this.w.close();
			return 0;
			} 
		catch (final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.w);
			CloserUtil.close(stream);
			CloserUtil.close(samReader);
			this.w=null;
			}
		}

	
	public static void main(final String[] args) {
		new Biostar172515().instanceMainWithExit(args);
	}

	}
