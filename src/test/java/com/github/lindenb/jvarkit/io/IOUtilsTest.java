package com.github.lindenb.jvarkit.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.zip.GZIPInputStream;

import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;
import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.lang.StringUtilsTest;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.iterator.LineIteratorTest;

@AlsoTest({LineIteratorTest.class,StringUtilsTest.class})
public class IOUtilsTest {
	private final TestSupport support = new TestSupport();
	
	private void write(OutputStream os) throws IOException {
		for(int i=0;i<10;i++) os.write((byte)'A');
	}
	
	private void read(InputStream is) throws IOException {
		for(int i=0;i<10;i++) Assert.assertEquals((int)'A', is.read());
		Assert.assertEquals(-1, is.read());
	}
	
	
@Test
public void testCompressed() {
	Assert.assertTrue(IOUtils.isCompressed(Paths.get("x.gz")));
	Assert.assertTrue(IOUtils.isCompressed(Paths.get("x.bgz")));
	Assert.assertTrue(IOUtils.isCompressed(Paths.get("x.bz2")));
	Assert.assertFalse(IOUtils.isCompressed(Paths.get("x.vcf")));
}
	
@Test
void testWrite() throws IOException {
	try {
		Path p = support.createTmpPath(".vcf.gz");
		OutputStream os = IOUtils.openPathForWriting(p);
		write(os);
		os.close();
		
		InputStream in = new GZIPInputStream(Files.newInputStream(p));
		read(in);
		in.close();
		
		in = IOUtils.openPathForReading(p);
		in.close();
		
		p = support.createTmpPath(".txt.gz");
		os = IOUtils.openPathForWriting(p);
		write(os);
		os.close();
		
		in = new GZIPInputStream(Files.newInputStream(p));
		in.close();
		in = IOUtils.openPathForReading(p);
		read(in);
		in.close();
		
		p = support.createTmpPath(".txt");
		os = IOUtils.openPathForWriting(p);
		write(os);
		os.close();
		
		in = Files.newInputStream(p);
		read(in);
		in.close();
		in = IOUtils.openPathForReading(p);
		read(in);
		in.close();
		os.close();
		
		p = support.createTmpPath(".bz2");
		os = IOUtils.openPathForWriting(p);
		write(os);
		os.close();
		
		in =  new BZip2CompressorInputStream(Files.newInputStream(p));
		read(in);
		in.close();
		in = IOUtils.openPathForReading(p);
		read(in);
		in.close();
		
		
	} finally {
		support.removeTmpFiles();
	}
	}


@Test
public void testRemoteInputStream() throws IOException {
	String url="https://en.wikipedia.org/wiki/Main_Page";
	Assert.assertTrue(IOUtils.isRemoteURI(url));
	boolean got_byte=false;
	try(InputStream br=IOUtils.openURIForReading(url))
		{
		;
		got_byte  = (br.read()!=-1);
		}
	Assert.assertTrue(got_byte);
	}


@Test
public void testRemoteReader() throws IOException {
	String url="https://en.wikipedia.org/wiki/Main_Page";
	Assert.assertTrue(IOUtils.isRemoteURI(url));
	boolean got_byte=false;
	try(BufferedReader br=IOUtils.openURIForBufferedReading(url))
		{
		got_byte    = (br.read()!=-1);
		}
	Assert.assertTrue(got_byte);
	}
@Test
public void testFileSuffix() {
	Assert.assertEquals(IOUtils.getFileSuffix(Paths.get("x.vcf")),".vcf");
	Assert.assertNotEquals(IOUtils.getFileSuffix(Paths.get("x.vcf")),"vcf");
	}
}
