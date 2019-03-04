package com.github.lindenb.jvarkit.io;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class ArchiveFactoryTest {
	private final TestSupport support = new TestSupport();

	private void fillArchive(ArchiveFactory archive) throws IOException {
		PrintWriter pw=archive.openWriter("jeter.txt");
		pw.print("ok");
		pw.flush();
		pw.close();
	
	}
	
@Test
public void testDirectory() throws IOException {
	Path dir = Files.createTempDirectory("tmp.");
	try {
		ArchiveFactory archive =  ArchiveFactory.open(dir);
		fillArchive(archive);
		archive.close();
		Assert.assertTrue(Files.exists(dir.resolve("jeter.txt")));
	} finally
	{
		Files.delete(dir.resolve("jeter.txt"));
		Files.delete(dir);
	}
}

@Test
public void testZip() throws IOException {
	try {
		Path zipFile = support.createTmpPath(".zip");
		ArchiveFactory archive =  ArchiveFactory.open(zipFile);
		fillArchive(archive);
		archive.close();
		support.assertZip(zipFile);
		}
	finally
		{
		support.removeTmpFiles();
		}
	}

}
