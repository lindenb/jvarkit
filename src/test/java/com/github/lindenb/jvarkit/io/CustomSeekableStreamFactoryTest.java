package com.github.lindenb.jvarkit.io;

import java.io.IOException;

import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;

public class CustomSeekableStreamFactoryTest  {
	private final TestSupport support = new TestSupport();

	private CustomSeekableStreamFactory getFactory() {
		return new CustomSeekableStreamFactory(SeekableStreamFactory.getInstance());
	}
	
	private void consumme(final SeekableStream st)throws IOException
		{
		st.seek(10);
		byte array[]=new byte[100];
		st.readFully(array);
		}
	
	@Test
	public void testFile() throws IOException {
		SeekableStream st = getFactory().getStreamFor(support.resource("S5.bam"));
		consumme(st);
		st.close();
		}
	@Test
	public void testURL() throws IOException {
		SeekableStream st = getFactory().getStreamFor("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/release/2009_02/Pilot1/md5.txt");
		consumme(st);
		st.close();
		}
	
}
