package com.github.lindenb.jvarkit.util.picard;

import java.io.File;

import org.testng.annotations.Test;

public class SAMSequenceDictionaryFactoryTest {
	@Test()
	public void loadExistingDict() throws Exception
		{
		load("Homo_sapiens_assembly18.trimmed.fasta");
		}
	@Test()
	public void loadNonExistingDict() throws Exception
		{
		//faidx exists, no dict
		load("Homo_sapiens_assembly18.trimmed.nodict.fasta");
		}
	
	private void load(String filename)throws Exception
		{
		new SAMSequenceDictionaryFactory().load(new File("htsjdk/testdata/htsjdk/samtools/reference/"+filename));
		}
}
