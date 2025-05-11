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

import java.io.IOException;

import com.github.lindenb.jvarkit.fastq.FastqPairedWriter;
import com.github.lindenb.jvarkit.fastq.FastqRecordPair;
import com.github.lindenb.jvarkit.jcommander.OnePassFastqLauncher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.CloseableIterator;


@Program(name="fastqphred64to33",
	keywords={"fastq"},
	description="Convert Illumina Fastq 64 encoding to Fastq 33"
	)
public class ConvertPhred64toFastq33 extends OnePassFastqLauncher
	{
	private static final Logger LOG = Logger.of(ConvertPhred64toFastq33.class);


	@Override
	protected Logger getLogger() {
		return LOG;
		}

	@Override
	protected int runPairedEnd(CloseableIterator<FastqRecordPair> iter, FastqPairedWriter fws) throws IOException {
		while(iter.hasNext()) {
			final FastqRecordPair pair = iter.next();
			fws.write(convert(pair.get(0)),convert(pair.get(1)));
			}
		return 0;
		}
	
	@Override
	protected int runSingleEnd(FastqReader iter, FastqWriter fws) throws IOException {
		while(iter.hasNext()) {
			fws.write(convert(iter.next()));
			}
		return 0;
		}
	
	private FastqRecord convert(FastqRecord rec) throws IOException
		{
		byte quals[]=rec.getBaseQualityString().getBytes();
		for(int i=0;i< quals.length;++i )
			{
			quals[i]=(byte)(quals[i]-64+33);
			if(quals[i]<33 || quals[i]>126)
				{
				throw new IOException("q="+(int)quals[i]);
				}
			}
		String name=rec.getReadName();
		int diez=name.indexOf('#');
		if(diez!=-1) name=name.substring(0, diez);
		
		return new FastqRecord(
				name,
				rec.getReadBases(),
				rec.getBaseQualityHeader() == null || rec.getReadName().equals(rec.getBaseQualityHeader())? "" : rec.getBaseQualityHeader(),
				quals
				);
					
		}
	

	public static void main(final String[] args) {
		new ConvertPhred64toFastq33().instanceMainWithExit(args);
	}

}
