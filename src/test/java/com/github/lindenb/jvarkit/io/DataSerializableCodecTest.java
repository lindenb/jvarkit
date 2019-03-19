package com.github.lindenb.jvarkit.io;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;


public class DataSerializableCodecTest {
	private static class XClass implements DataSerializable
		{
		String s;
		int i;
		@Override
		public void readDataInputStream(DataInputStream in) throws IOException {
			s=in.readUTF();
			i=in.readInt();
		}

		@Override
		public void writeDataOutputStream(DataOutputStream daos) throws IOException {
			daos.writeUTF(s);
			daos.writeInt(i);
		}
		
		}
@Test
void test01() throws IOException{
	ByteArrayOutputStream baos=new ByteArrayOutputStream();
	DataOutputStream daos = new DataOutputStream(baos);
	for(int i=0;i< 100;i++) {
		XClass x=new XClass();
		x.s=""+i;
		x.i=i;
		x.writeDataOutputStream(daos);
		}
	daos.flush();
	baos.flush();
	ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
	DataInputStream dais=new DataInputStream(bais);
	for(int i=0;i< 100;i++) {
		XClass x=new XClass();
		x.readDataInputStream(dais);
		Assert.assertEquals(x.s,""+i);
		Assert.assertEquals(x.i,i);
		}
	}

@Test
void test02() throws IOException{
	
	DataSerializableCodec<XClass> codec=new DataSerializableCodec<>(()->new XClass());
	ByteArrayOutputStream baos=new ByteArrayOutputStream();
	DataOutputStream daos = new DataOutputStream(baos);
	for(int i=0;i< 100;i++) {
		XClass x=new XClass();
		x.s=""+i;
		x.i=i;
		codec.encode(daos, x);
		}
	daos.flush();
	baos.flush();

	
	
	ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
	DataInputStream dais=new DataInputStream(bais);
	for(int i=0;i< 100;i++) {
		XClass x=codec.decode(dais);
		Assert.assertNotNull(x);
		Assert.assertEquals(x.s,""+i);
		Assert.assertEquals(x.i,i);
		}
	}


}
