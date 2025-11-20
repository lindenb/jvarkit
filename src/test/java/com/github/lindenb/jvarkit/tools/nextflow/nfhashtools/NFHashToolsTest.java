package com.github.lindenb.jvarkit.tools.nextflow.nfhashtools;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class NFHashToolsTest {

private Path createJSON(TestSupport support,String value) throws IOException {
    final Path logfile = support.createTmpPath(".log");
    Files.writeString(logfile, 
    		"Nov-19 22:23:53.071 [Actor Thread 96] INFO  nextflow.processor.TaskProcessor - [TASK1] cache hash: 0d4bb965b6140f6fc8a87f840c295e30; mode: LENIENT; entries: [\n" +
    		"    {" + 
    		"        \"hash\": \"5350de1bf7fd1f286a1ad5baa7af30f5\"," + 
    		"        \"type\": \"java.util.UUID\", " +
    		"        \"value\": \""+value+"\" "  +
    		"    }\n"+
    		"]" 
    		
    		);
    return logfile;
    
	}
	
@Test
public void testOneFile()throws IOException {
	final TestSupport support= new TestSupport();
	try {
		Path p1 = createJSON(support,"V1");
		Path out = support.createTmpPath(".json");
		Assert.assertEquals(new NFHashTools().instanceMain(new String[] {
				"-o",out.toString(),
				p1.toString()
				}),0);
		Assert.assertEquals(new NFHashTools().instanceMain(new String[] {
				"-o",out.toString(),
				"-t","V1",
				p1.toString()
				}),0);
		}
	finally
		{
		support.removeTmpFiles();
		}
	}
@Test
public void testTwoFiles()throws IOException {
	final TestSupport support= new TestSupport();
	try {
		Path p1 = createJSON(support,"V1");
		Path p2 = createJSON(support,"V2");
		Path out = support.createTmpPath(".json");
		Assert.assertEquals(new NFHashTools().instanceMain(new String[] {
				"-o",out.toString(),
				p1.toString(),
				p2.toString()
				}),0);
		
		}
	finally
		{
		support.removeTmpFiles();
		}
	}
}
