package com.github.lindenb.jvarkit.tools.genbank;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.ncbi.NcbiApiKey;
import com.github.lindenb.jvarkit.util.ncbi.NcbiConstants;

public class GenBankToGff3Test{
	private final TestSupport support = new TestSupport();
	private final NcbiApiKey ncbiApiKey = new NcbiApiKey();
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
			{"AF338247.1,D38149.1"},//rota + rota
			{"KU956010.1"},//rota
			{"NC_001604.1"},// Enterobacteria phage T7, complete genome.
			{"NC_014415.1"}//small eukaryote
			};
	}
	
	private void testXml(final Path gbxml) throws IOException
		{
		final Path out = support.createTmpPath(".txt");
		Assert.assertEquals(
			new GenbankToGff3().instanceMain(new String[] {
			"-o",
			out.toString(),
			gbxml.toString()
			}),0);
		support.assertIsNotEmpty(out);
		}
	
	@Test(dataProvider="src1")
	public void test1(final String acns) throws IOException {
		try {
			final Path gbxml = support.createTmpPath(".xml");
			final InputStream xmlin = IOUtils.openURIForReading(NcbiConstants.efetch()+
					"?db=nuccore&id="+acns+"&rettype=gbwithparts&retmode=xml"
							+ this.ncbiApiKey.getAmpParamValue());
					
			IOUtils.copyTo(xmlin, gbxml);
			xmlin.close();
			support.assertIsXml(gbxml);
			testXml(gbxml);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}
	
	@DataProvider(name = "src2")
	public Object[][] createData2() {
		return new Object[][]{
			{support.resource("rotavirus_rf.gb.xml")}
			};
	}
	
	@Test(dataProvider="src2")
	public void test2(final String xmlfile) throws IOException {
		try {
			testXml(Paths.get(xmlfile));
			}
		finally
			{
			support.removeTmpFiles();
			}
		}
	
	
}
