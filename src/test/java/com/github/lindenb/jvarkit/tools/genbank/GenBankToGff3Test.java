package com.github.lindenb.jvarkit.tools.genbank;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.TestUtils;
import com.github.lindenb.jvarkit.util.ncbi.NcbiConstants;

public class GenBankToGff3Test extends TestUtils{
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
			{"AF338247.1,D38149.1"},//rota + rota
			{"KU956010.1"},//rota
			{"NC_001604.1"},// Enterobacteria phage T7, complete genome.
			{"NC_014415.1"}//small eukaryote
			};
	}
		
	@Test(dataProvider="src1")
	public void test1(final String acns) throws IOException {
		final File gbxml = createTmpFile(".xml");
		final InputStream xmlin = IOUtils.openURIForReading(NcbiConstants.efetch()+
				"?db=nuccore&id="+acns+"&rettype=gbwithparts&retmode=xml"
						+ super.ncbiApiKey.getAmpParamValue());
				
		IOUtils.copyTo(xmlin, gbxml);
		xmlin.close();
		assertIsXml(gbxml);
		
		final File out = createTmpFile(".txt");
		Assert.assertEquals(
			new GenbankToGff3().instanceMain(newCmd().
			add("-o").add(out).
			add(gbxml).
			make()
			),0);
		assertIsNotEmpty(out);
		}
}
