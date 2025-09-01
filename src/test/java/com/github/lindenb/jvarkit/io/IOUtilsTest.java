package com.github.lindenb.jvarkit.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.zip.GZIPInputStream;

import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;
import org.testng.Assert;
import org.testng.annotations.Test;



import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class IOUtilsTest {
	private static final String STABLE_URL="https://mastodon.social//api/v1/timelines/public";
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
		try(OutputStream os = IOUtils.openPathForWriting(p)) {
		write(os);
		}
		
		try(InputStream in = new GZIPInputStream(Files.newInputStream(p))) {
		read(in);
		}
		
		try(InputStream in = IOUtils.openPathForReading(p)) {
		}
		
		p = support.createTmpPath(".txt.gz");
		try(OutputStream os = IOUtils.openPathForWriting(p)) {
		write(os);
		}
		
		try(InputStream in = new GZIPInputStream(Files.newInputStream(p))) {
		}
		
		
		try(InputStream in = IOUtils.openPathForReading(p)) {
		read(in);
		}
		
		p = support.createTmpPath(".txt");
		try(OutputStream os = IOUtils.openPathForWriting(p)) {
		write(os);
		}
		
		try(InputStream in  = Files.newInputStream(p)) {
		read(in);
		}
		
		try(InputStream in = IOUtils.openPathForReading(p)) {
		read(in);
		}
		
		p = support.createTmpPath(".bz2");
		try(OutputStream os = IOUtils.openPathForWriting(p)) {
			write(os);
			}
		
		try(InputStream in =  new BZip2CompressorInputStream(Files.newInputStream(p))) {
			read(in);
			}

		try(InputStream in =IOUtils.openPathForReading(p)) {
			read(in);
			}
			
		
	} finally {
		support.removeTmpFiles();
	}
	}


@Test
public void testRemoteInputStream() throws IOException {
	String url=STABLE_URL;
	Assert.assertTrue(IOUtils.isRemoteURI(url));
	boolean got_byte=false;
	try(InputStream br=IOUtils.openURIForReading(url))
		{
		got_byte  = (br.read()!=-1);
		}
	Assert.assertTrue(got_byte);
	}


@Test
public void testRemoteReader() throws IOException {
	String url=STABLE_URL;
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

@Test
public void testURL() throws Exception {
	String str=STABLE_URL;
	Assert.assertEquals(IOUtils.toURL(str),new URI(str).toURL());
	}

}
