package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMFileHeader.SortOrder;

public class Biostar78400Test {
	private final TestSupport support = new TestSupport();

	
@Test
public void test01() throws IOException {
	try {
	final String flowcell="HS20001259127";
	final String lane = "1";
	final Path in = support.createTmpPath(".bam");
	SAMFileHeader header=new SAMFileHeader();
	header.setSortOrder(SortOrder.unsorted);
	SAMFileWriter sfw =new SAMFileWriterFactory().makeBAMWriter(header, true, in);
	DefaultSAMRecordFactory recfactory = new DefaultSAMRecordFactory();
	SAMRecord rec=recfactory.createSAMRecord(header);
	rec.setReadName(flowcell+":"+lane+":1210:15640:52255");
	rec.setReadString("GAATTC");
	rec.setBaseQualityString("222222");
	SAMUtils.makeReadUnmapped(rec);
	sfw.addAlignment(rec);
	sfw.close();
	support.assertIsValidBam(in);
	
	final Path xml = support.createTmpPath(".xml");
	PrintWriter pw = new PrintWriter(Files.newOutputStream(xml));
	pw.println("<?xml version=\"1.0\"?><read-groups>"
			+ "<flowcell name=\""+flowcell+"\"><lane index=\""+ lane+"\">"
			+ "<group ID=\"X1\"><library>L1</library><platform>P1</platform>"
			+ "<sample>S1</sample><platformunit>PU1</platformunit>"
			+ "<center>C1</center><description>blabla</description></group>"
			+ "</lane></flowcell><flowcell name=\"HS20001259128\">"
			+ "<lane index=\"2\"><group ID=\"x2\"><library>L2</library>"
			+ "<platform>P2</platform><sample>S2</sample><platformunit>PU1</platformunit>"
			+ "<center>C1</center><description>blabla</description></group></lane>"
			+ "</flowcell></read-groups>"
			);
	pw.flush();
	pw.close();
	support.assertIsXml(xml);
		
		final Path out = support.createTmpPath(".bam");
		Assert.assertEquals(
			new Biostar78400().instanceMain(new String[] {
			"-o",out.toString(),
			"-x",xml.toString(),
			in.toString()
			}),0);
	
		SamReader r= SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT).open(out);
		Assert.assertTrue(r.getFileHeader()!=null);
		Assert.assertTrue(r.getFileHeader().getReadGroups()!=null);
		Assert.assertFalse(r.getFileHeader().getReadGroups().isEmpty());
		SAMRecordIterator iter = r.iterator();
		Assert.assertTrue(iter.hasNext());
		rec = iter.next();
		SAMReadGroupRecord rg = rec.getReadGroup();
		Assert.assertNotNull(rg);
		Assert.assertEquals(rg.getId(),"X1");
		Assert.assertEquals(rg.getSample(),"S1");
		
		Assert.assertFalse(iter.hasNext());
		
		iter.close();
		r.close();
		} 
	finally {
		support.removeTmpFiles();
		}
	}
	
}
