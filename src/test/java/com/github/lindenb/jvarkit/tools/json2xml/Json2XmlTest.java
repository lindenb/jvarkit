package com.github.lindenb.jvarkit.tools.json2xml;

import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class Json2XmlTest {
	private Path writeRandomJSON(TestSupport support) throws Exception {
		Path p = support.createTmpPath(".json");
		Files.writeString(p, "");
		return p;
		}
	@Test
	public void withoutNS() throws Exception {
		final TestSupport support= new TestSupport();
		try {
			Path p = writeRandomJSON(support);
			Path output = support.createTmpPath(".xml");
			Assert.assertEquals(new Json2Xml().instanceMain(new String[]{
	        		"-o",output.toString(),
	        		"--ns","",
	        		p.toString()
	        	}),0);
	        support.assertIsXml(output);
		}
		
		finally {
			support.removeTmpFiles();
		}
	}
	
	@Test
	public void withNS() throws Exception {
		final TestSupport support= new TestSupport();
		try {
			Path p = writeRandomJSON(support);
			Path output = support.createTmpPath(".xml");
			Assert.assertEquals(new Json2Xml().instanceMain(new String[]{
	        		"-o",output.toString(),
	        		p.toString()
	        	}),0);
	        support.assertIsXml(output);
		}
		
		finally {
			support.removeTmpFiles();
		}
	}

}
