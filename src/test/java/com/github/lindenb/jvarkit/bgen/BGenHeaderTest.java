package com.github.lindenb.jvarkit.bgen;


import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class BGenHeaderTest extends TestSupport{
	private BGenHeader readHeader(String f) throws IOException {
		BGenCodec codec = new BGenCodec();
		final Path p =Paths.get(SRC_TEST_RESOURCE+File.separator+"BGEN"+File.separator+f+BGenUtils.FILE_SUFFIX);
		Assert.assertTrue(Files.exists(p));
		Assert.assertTrue(Files.isReadable(p));
		try(BGenUtils.RandomAccessStream in=codec.makeSourceFromStream(Files.newInputStream(p))) {
			return codec.readActualHeader(in);
			}
		}
	
	private void checkSample(BGenHeader hdr,int expect) {
		Assert.assertNotNull(hdr.getSamples());
		Assert.assertEquals(hdr.getNSamples(),expect);
		for(int i=0;i< hdr.getNSamples();++i) {
			String sn = hdr.getSamples().get(i);
			Assert.assertTrue(sn.trim().length()>0);
			Assert.assertEquals(i,hdr.getSampleIndex(sn));
			}
		}
	
	@Test
	public void testReadHeaderV11() throws IOException {
		BGenHeader hdr = readHeader("example.001.no_sample_block.1.1");
		Assert.assertNotNull(hdr);
		Assert.assertEquals(hdr.getLayout(),BGenUtils.Layout.LAYOUT_1);
		Assert.assertEquals(hdr.getCompression(),BGenUtils.Compression.ZLIB);
		Assert.assertTrue(hdr.hasAnonymousSamples());
		checkSample(hdr,11);
		}
}
