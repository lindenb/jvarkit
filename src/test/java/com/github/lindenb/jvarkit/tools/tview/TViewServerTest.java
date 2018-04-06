package com.github.lindenb.jvarkit.tools.tview;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URLEncoder;
import java.util.Arrays;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.util.CloserUtil;

public class TViewServerTest  extends TestUtils {
	
	@DataProvider(name="rf_regions")
	public Object[][] getDataRegions() throws IOException  {
		return randomIntervalsFromDict(new File(SRC_TEST_RESOURCE+"/rotavirus_rf.dict"),3).
				stream().
				map(I->I.getContig()+":"+I.getStart()+"-"+I.getEnd()).
				map(S->new Object[]{S}).
				toArray(x->new Object[x][])
				;
		}
		
	@Test(dataProvider="rf_regions")
	public void test01(final String rgn) throws IOException {
		final int port = 9090;
		final File htmlOut = super.createTmpFile(".html");
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
		
		new  TViewServer().instanceMain(newCmd().add(
				"-P",port,
				"-R",SRC_TEST_RESOURCE+"/rotavirus_rf.fa",
				"--shutdown-after","60",
				SRC_TEST_RESOURCE+"/S1.bam",
				SRC_TEST_RESOURCE+"/S2.bam"
				).make());
		
		assertIsXml(htmlOut);
		}


}
