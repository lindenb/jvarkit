package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class Biostar86480Test {
	
	private final TestSupport support = new TestSupport();

	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		List<String> L1=  support.allFasta().collect(Collectors.toList());
		List<String> L2= Arrays.asList("-E HhaI","-E AluI","-E HhaI -E AluI");
		List<Object[]> L3=new ArrayList<>();
		for(String s1:L1) {
			for(String s2:L2) {
				L3.add(new Object[] {s1,s2});
			}
		}
		return (Object[][])L1.toArray();
		}
		
	@Test(dataProvider="src1")
	public void test1(final String infFasta,final String enz) throws IOException {
		try {
		final Path out = support.createTmpPath(".txt");
		final Biostar86480 cmd =new Biostar86480();
		List<String> args = new ArrayList<>();
		args.add("-o");
		args.add(out.toString());
		args.addAll(Arrays.asList(enz.split("[ \t]+")));
		args.add(infFasta);
		
		Assert.assertEquals(cmd.instanceMain(args),0);
		support.assertIsNotEmpty(out);
		} finally {
			support.removeTmpFiles();
		}
		}
}
