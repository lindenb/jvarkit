/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.dict;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.Writer;
import java.math.BigInteger;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceDictionaryCodec;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.util.ParsingUtils;

public class DictionaryCreator {

public static void createFromFasta(Path fasta) throws IOException {
	final Path dictf = ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(fasta);
	final SAMSequenceDictionary dict= new DictionaryCreator().fromFasta(fasta);
	try(Writer w=Files.newBufferedWriter(dictf)) {
        new SAMSequenceDictionaryCodec(w).encode(dict);
        w.flush();
    	}
	}
	
public SAMSequenceDictionary fromFasta(URL url) throws IOException {
	try(InputStream in=IOUtils.mayBeGzippedInputStream(ParsingUtils.getURLHelper(url).openInputStream())) {
		return fromFasta(in);
		}
	}
	
	
public SAMSequenceDictionary fromFasta(Path path) throws IOException {
	try(InputStream in=IOUtil.openFileForReading(path)) {
		return fromFasta(in);
		}
	}

public SAMSequenceDictionary fromFasta(InputStream is) throws IOException {
	final List<SAMSequenceRecord> L=new ArrayList<>();
	
	try(BufferedInputStream in = new BufferedInputStream(is)) {
		String name=null;
		final MessageDigest md5;
		
		try {
            md5 = MessageDigest.getInstance("MD5");
            md5.reset();
       		}
        catch (NoSuchAlgorithmException e) {
            throw new RuntimeException("MD5 algorithm not found", e);
        	}
		
		int length=0;
		boolean at_start=true;
		for(;;)
				{
				int c=in.read();
		        if (c==-1 || (at_start && c== '>')) {
		        	if(name!=null)
		        		{
		        		if(length==0)throw new IOException("empty sequence for "+name);
		        		final SAMSequenceRecord ssr=new SAMSequenceRecord(name, length);
		        		

		        		String hash = new BigInteger(1, md5.digest()).toString(16);
		                while(hash.length() < 32) {
		                	hash+='0';
		                	}
		                ssr.setMd5(hash);
		        		
		        		L.add(ssr);
		        		}
		        	if(c==-1) break;
		        	length=0;
		        	name=null;
		        	final StringBuilder sb=new StringBuilder();
		        	while((c=in.read())!=-1 && c!='\n' && c!='\r') {
		        		sb.append((char)c);
		        		}
		        	if(c=='\n' || c=='\r') {
		            	at_start=true;
		            	}
		        	name=SAMSequenceRecord.truncateSequenceName(sb.toString());
		        	if(StringUtils.isBlank(name)) throw new IOException("empty sequence name");
		        	md5.reset();
		        	}
		        else if(c=='\n' || c=='\r') {
		        	at_start=true;
		        	}
		        else if(name!=null)
		        	{
		        	md5.update((byte)c);
		        	length++;
		        	}
		        else if(!Character.isWhitespace(c))
		        	{
		        	throw new IOException("Illegal character "+(char)c);
		        	}
		        }/* end for */
	        }
	return new SAMSequenceDictionary(L);
	}
	
}
