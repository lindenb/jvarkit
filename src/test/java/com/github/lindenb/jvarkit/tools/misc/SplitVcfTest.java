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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class SplitVcfTest extends TestUtils {

@Test
public void testWithFile() throws IOException {
	File in = super.createTmpFile(".txt");
	final PrintWriter pw = new PrintWriter(in);
	pw.println("X1\tRF01:970-970");
	pw.println("X2\tRF02:578-578 RF02:1962-1962");
	pw.flush();
	pw.close();
	
	final String template = IOUtils.getDefaultTmpDir().getPath()+File.separator+ "__GROUPID__.vcf.gz";
	
	Assert.assertEquals(new SplitVcf().instanceMain(new String[] {
			"--unmapped","UNMAPPED",
			"-g",in.getPath(),
			"-o",template,
			SRC_TEST_RESOURCE+"/rotavirus_rf.vcf.gz"
			}),0);
	in = new File(template.replaceAll("__GROUPID__", "X1"));
	Assert.assertEquals(super.variantStream(in).count(), 1L);
	Assert.assertTrue(in.delete());
	
	in = new File(template.replaceAll("__GROUPID__", "X2"));
	Assert.assertEquals(super.variantStream(in).count(), 2L);
	Assert.assertTrue(in.delete());
	
	in = new File(template.replaceAll("__GROUPID__", "UNMAPPED"));
	Assert.assertTrue(super.variantStream(in).count()>0);
	Assert.assertTrue(in.delete());
	}
@Test
public void testWithDict() throws IOException {
	final String template = IOUtils.getDefaultTmpDir().getPath()+File.separator+ "__GROUPID__.vcf.gz";
	final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(new File(SRC_TEST_RESOURCE+"/rotavirus_rf.dict"));
	Assert.assertEquals(new SplitVcf().instanceMain(new String[] {
			"--unmapped","UNMAPPED",
			"-o",template,
			SRC_TEST_RESOURCE+"/rotavirus_rf.vcf.gz"
			}),0);
	File in;
	for(final SAMSequenceRecord rec:dict.getSequences()) {
		in = new File(template.replaceAll("__GROUPID__",rec.getSequenceName()));
		super.assertIsVcf(in);
		Assert.assertTrue(in.delete());
		}
	
	in = new File(template.replaceAll("__GROUPID__", "UNMAPPED"));
	if(in.exists()) in.delete();
	}

}
