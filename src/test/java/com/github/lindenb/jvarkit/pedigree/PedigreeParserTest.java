package com.github.lindenb.jvarkit.pedigree;

import java.io.IOException;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;
/*
FAM01	S1	0	0	1	1
FAM01	S2	0	0	2	1
FAM01	S3	S1	S2	1	1
FAM02	S4	0	0	1	0
FAM02	S5	0	0	2	0
FAM02	S6	S4	S5	2	0

 */
public class PedigreeParserTest {
	
private final TestSupport support = new TestSupport();
	
	
@Test
public void test01() throws IOException{
	final PedigreeParser parser = new PedigreeParser();
	final Pedigree ped=parser.parse(Paths.get(support.resource("test_vcf01.ped")));
	Assert.assertNotNull(ped);
	Assert.assertFalse(ped.isEmpty());
	Assert.assertEquals(ped.getSamples().size(),6);
	Assert.assertEquals(ped.getTrios().size(),2);
	Assert.assertEquals(ped.getFamilies().size(),2);
	Assert.assertNull(ped.getFamilyById("xx"));
	Assert.assertEquals(ped.getTrios().size(),2);

	
	Family fam = ped.getFamilyById("FAM01");
	Assert.assertNotNull(fam);
	
	
	Assert.assertEquals(fam.getSamples().size(),3);
	Sample sample=fam.getSampleById("xxx");
	Assert.assertNull(sample);
	sample=fam.getSampleById("S3");
	Assert.assertNotNull(sample);
	
	Assert.assertEquals(fam.getTrios().size(),1);
	Assert.assertEquals(fam.getTrios().iterator().next().getChild(),sample);
	
	
	Assert.assertEquals(sample.getSex(),Sex.male);
	Assert.assertTrue(sample.isMale());
	Assert.assertFalse(sample.isFemale());

	Assert.assertEquals(sample.getStatus(),Status.affected);
	Assert.assertTrue(sample.isAffected());
	Assert.assertFalse(sample.isUnaffected());
	
	Assert.assertTrue(sample.hasFather());
	Assert.assertNotNull(sample.getFather());
	Assert.assertEquals(sample.getFather().getId(),"S1");
	Assert.assertTrue(sample.getFather().getParents().isEmpty());

	
	Assert.assertTrue(sample.hasMother());
	Assert.assertNotNull(sample.getMother());
	Assert.assertEquals(sample.getMother().getId(),"S2");
	Assert.assertTrue(sample.getMother().getParents().isEmpty());
	}
}
