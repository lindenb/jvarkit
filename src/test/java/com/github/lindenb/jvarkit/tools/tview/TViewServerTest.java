package com.github.lindenb.jvarkit.tools.tview;
import java.io.IOException;
import java.io.InputStream;
import java.net.URLEncoder;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

import htsjdk.samtools.util.CloserUtil;

@AlsoTest(LauncherTest.class)
public class TViewServerTest {
	
private final TestSupport support = new TestSupport();

	
	@DataProvider(name="rf_regions")
	public Object[][] getDataRegions() throws IOException  {
		return support.toArrayArray(
				support.randomIntervalsFromDict(Paths.get(support.resource("rotavirus_rf.dict")),10,1000).
				stream().
				map(I->I.getContig()+":"+I.getStart()+"-"+I.getEnd()).
				map(S->new Object[]{S})
				)
				;
		}
		
	@Test(dataProvider="rf_regions")
	public void test01(final String rgn) throws IOException {
		try {
			final int port = 9090;
			final Path htmlOut = support.createTmpPath(".html");
			final String uri = "http://localhost:"+port+"/?rgn="+URLEncoder.encode(rgn, "UTF-8");
			new java.util.Timer().schedule( 
			        new java.util.TimerTask() {
			            @Override
			            public void run() {
			            	InputStream in= null;
			            	try {
			            		 in=IOUtils.openURIForReading(uri);
			            		IOUtils.copyTo(in, htmlOut);
			            	} catch(IOException err) {
			            		Assert.fail("cannot read "+uri, err);
			            	}
			            	finally {
			            		CloserUtil.close(in);
			            	}
			            }
			        },30);
			
			new  TViewServer().instanceMain(new String[] {
					"-P",String.valueOf(port),
					"-R",support.resource("rotavirus_rf.fa"),
					"--shutdown-after","60",
					support.resource("S1.bam"),
					support.resource("S2.bam")
				});
			
			support.assertIsXml(htmlOut);
		} finally {
			support.removeTmpFiles();
			}
		}


}
