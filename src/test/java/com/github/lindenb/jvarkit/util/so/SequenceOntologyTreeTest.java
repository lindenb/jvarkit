/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


*/
package com.github.lindenb.jvarkit.util.so;

import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URL;

import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


public class SequenceOntologyTreeTest
	{
	private SequenceOntologyTree owlTree;
	@BeforeClass
	public void downloadSoOwl() throws IOException {
		try {
		URL url=new URL("https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/so-simple.owl");
		HttpURLConnection con = (HttpURLConnection)url.openConnection();
		HttpURLConnection.setFollowRedirects(true);
		con.connect();
		InputStream in = con.getInputStream();
		this.owlTree = SequenceOntologyTree.createFromInputStream(in);
		in.close();
		con.disconnect();
		}
		catch(final Throwable err) {
			owlTree = SequenceOntologyTree.createDefault();
		}
	}
	@AfterClass
	public void disposeOwl()
		{
		owlTree=null;
		}
	
	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
	 return new Object[][] {
	   { "SO:0001583","missense_variant","SO:0001818" },
	   { "SO:0001818","protein_altering_variant","SO:0001576" }
	 };
	}

	
	 
	@Test(dataProvider="src1")
	public void test1(final String acn,final String label,final String parentAcn) 
		{
		testTerms(SequenceOntologyTree.getInstance(),acn,label,parentAcn);
		}
	
	@Test(dataProvider="src1")
	public void test2(final String acn,final String label,final String parentAcn) 
		throws IOException
		{
		testTerms(this.owlTree,acn,label,parentAcn);
		}
	
	
	
	private void testTerms(final SequenceOntologyTree tree,final String acn,final String label,final String parentAcn) 
		{
		SequenceOntologyTree.Term t1 = tree.getTermByAcn(acn);
		Assert.assertNotNull(t1);
		Assert.assertEquals(t1.getAcn(), acn);
		Assert.assertEquals(t1.getLabel(), label);
		SequenceOntologyTree.Term t2 = tree.getTermByLabel(label);
		Assert.assertEquals(t1,t2);
		
		SequenceOntologyTree.Term t3 = tree.getTermByAcn(parentAcn);
		Assert.assertNotNull(t3);
		Assert.assertEquals(t3.getAcn(), parentAcn);

		
		Assert.assertTrue(t3.getAllDescendants().contains(t1));
		Assert.assertFalse(t1.getAllDescendants().contains(t3));
		
		Assert.assertTrue(t2.isChildrenOf(t3));
		Assert.assertFalse(t3.isChildrenOf(t2));

		}
	 

	}
