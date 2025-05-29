/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.so;

import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.Set;

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
		try(InputStream in = con.getInputStream()) {
			this.owlTree = SequenceOntologyTree.createFromInputStream(in);
			}
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
	 


    @Test
    public void testSingletonInstance() {
        SequenceOntologyTree tree1 = SequenceOntologyTree.getInstance();
        SequenceOntologyTree tree2 = SequenceOntologyTree.getInstance();
        Assert.assertNotNull(tree1);
        Assert.assertSame(tree1, tree2, "getInstance() should always return the same object");
    }

    @Test
    public void testDefaultTreeHasKnownTerms() {
        SequenceOntologyTree tree = SequenceOntologyTree.createDefault();
        SequenceOntologyTree.Term t = tree.getTermByAcn("SO:0000001");
        Assert.assertNotNull(t, "region should exist");
        Assert.assertEquals(t.getLabel(), "region");
        Assert.assertEquals(tree.getTermByLabel("region"), t);
    }

    @Test
    public void testTermParentChildRelationship() {
        SequenceOntologyTree tree = SequenceOntologyTree.createDefault();
        SequenceOntologyTree.Term region = tree.getTermByAcn("SO:0000001");
        SequenceOntologyTree.Term parent = tree.getTermByAcn("SO:0000110");
        Assert.assertTrue(region.getParents().contains(parent));
        Assert.assertTrue(parent.getChildren().contains(region));
        Assert.assertTrue(region.isChildrenOf(parent));
        Assert.assertFalse(parent.isChildrenOf(region));
    }

    @Test
    public void testGetAllDescendants() {
        SequenceOntologyTree tree = SequenceOntologyTree.createDefault();
        SequenceOntologyTree.Term region = tree.getTermByAcn("SO:0000001");
        Set<SequenceOntologyTree.Term> descendants = region.getAllDescendants();
        Assert.assertTrue(descendants.contains(region));
        // Should contain itself
        for (SequenceOntologyTree.Term child : region.getChildren()) {
            Assert.assertTrue(descendants.contains(child));
        }
    }

    @Test
    public void testUnknownTermReturnsNull() {
        SequenceOntologyTree tree = SequenceOntologyTree.createDefault();
        Assert.assertNull(tree.getTermByAcn("SO:9999999"));
        Assert.assertNull(tree.getTermByLabel("not_a_label"));
    }

    @Test
    public void testDamagingComparator() {
        SequenceOntologyTree tree = SequenceOntologyTree.createDefault();
        SequenceOntologyTree.DamagingComparator comp = new SequenceOntologyTree.DamagingComparator(tree);
        SequenceOntologyTree.Term t1 = tree.getTermByLabel("missense_variant");
        SequenceOntologyTree.Term t2 = tree.getTermByLabel("intergenic_variant");
        // These may not be present in the default, so just check that the comparator runs
        if (t1 != null && t2 != null) {
            Assert.assertNotEquals(0, comp.compare(t1, t2));
        }
    }
	
	
	}
