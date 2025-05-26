package com.github.lindenb.jvarkit.liftover;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StringReader;

import org.testng.Assert;
import org.testng.annotations.Test;

public class LiftOverChainInputStreamTest {

	private String[]  mockChain() throws IOException {
	return	new String[]{"chain 1000 RF01 3302 + 0 2600 RF02 2687 + 0 2600 1","2600\t0\t0","0"};
		}
	
	@Test
	public void testConvertIdentity() throws IOException {
		String[] lines = mockChain();
		try(Reader in =new StringReader(String.join("\n",lines)+"\n")) {
			InputStream in2 =new LiftOverChainInputStream(in, null, null);
			try(BufferedReader br=new BufferedReader(new InputStreamReader(in2))) {
				for(int i=0;i< lines.length;i++) {
					String line=br.readLine();
					Assert.assertNotNull(line);
					Assert.assertEquals(line,lines[i]);
					}
				Assert.assertNull(br.readLine());
			}
		}
	}
	
	@Test
	public void testConvertChangeSrc() throws IOException {
		String[] lines = mockChain();
		try(Reader in =new StringReader(String.join("\n",lines)+"\n")) {
			InputStream in2 =new LiftOverChainInputStream(in, A->"chr"+A, null);
			try(BufferedReader br=new BufferedReader(new InputStreamReader(in2))) {
				for(int i=0;i< lines.length;i++) {
					String line=br.readLine();
					Assert.assertNotNull(line);
					if(i==0) {
						Assert.assertEquals(line,lines[i].replace("chain 1000 RF01", "chain 1000 chrRF01"));
						}
					else
						{
						Assert.assertEquals(line,lines[i]);
						}
					}
				Assert.assertNull(br.readLine());
			}
		}
	}
	
	@Test
	public void testConvertChangeDest() throws IOException {
		String[] lines = mockChain();
		try(Reader in =new StringReader(String.join("\n",lines)+"\n")) {
			InputStream in2 =new LiftOverChainInputStream(in,null,  A->"chr"+A);
			try(BufferedReader br=new BufferedReader(new InputStreamReader(in2))) {
				for(int i=0;i< lines.length;i++) {
					String line=br.readLine();
					Assert.assertNotNull(line);
					if(i==0) {
						Assert.assertEquals(line,lines[i].replace("2600 RF02", "2600 chrRF02"));
						}
					else
						{
						Assert.assertEquals(line,lines[i]);
						}
					}
				Assert.assertNull(br.readLine());
			}
		}
	}
		
		
		@Test
		public void testConvertNoMatch() throws IOException {
			String[] lines = mockChain();
			try(Reader in =new StringReader(String.join("\n",lines)+"\n")) {
				InputStream in2 =new LiftOverChainInputStream(in,null,  A->null);
				try(BufferedReader br=new BufferedReader(new InputStreamReader(in2))) {
					String line=br.readLine();
					Assert.assertNull(line);
				}
			}
		}
}
