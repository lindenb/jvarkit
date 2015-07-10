/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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


History:
* 2015 creation

*/
package com.github.lindenb.jvarkit.util.picard;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;

import java.io.File;
import java.io.OutputStream;

public enum SamFormat {
sam()
	{
	@Override
	public boolean isBinary() {
		return false;
		}
	@Override
		public SAMFileWriter openSAMFileWriterToFile(SAMFileWriterFactory sfw,
				File file, SAMFileHeader header, boolean presorted) {
		return sfw.makeSAMWriter(header, presorted, file);
		}
	@Override
		public SAMFileWriter openSAMFileWriterToStream(
				SAMFileWriterFactory sfw, OutputStream out,
				SAMFileHeader header, boolean presorted) {
		return sfw.makeSAMWriter(header, presorted, out);
		}
	},
bam()
	{
	@Override
	public boolean isBinary() {
		return true;
		}
	@Override
		public SAMFileWriter openSAMFileWriterToFile(SAMFileWriterFactory sfw,
				File file, SAMFileHeader header, boolean presorted) {
		return sfw.makeBAMWriter(header, presorted, file);
		}
	@Override
		public SAMFileWriter openSAMFileWriterToStream(
				SAMFileWriterFactory sfw, OutputStream out,
				SAMFileHeader header, boolean presorted) {
		return sfw.makeBAMWriter(header, presorted, out);
		}
	};
/** return file extension, prefixed with dot '.' */
public String getExtension()
	{
	return "."+this.name();
	}
/** tries to find the SAM format from file, returns null if not found */
public static SamFormat getSamFormatFromFile(File f)
	{
	return SamFormat.getSamFormatFromFile(f.getName());
	}	
/** tries to find the SAM format from file, returns null if not found */
public static SamFormat getSamFormatFromFile(String s)
	{
	for(SamFormat f:values())
		{
		if(s.endsWith(f.getExtension())) return f;
		}
	return null;
	}	
public  SAMFileWriter openSAMFileWriterToStdout(SAMFileWriterFactory sfw,SAMFileHeader header,boolean presorted)
	{
	return openSAMFileWriterToStream(sfw, System.out, header, presorted);
	}

public abstract SAMFileWriter openSAMFileWriterToStream(SAMFileWriterFactory sfw,OutputStream out,SAMFileHeader header,boolean presorted);
public abstract SAMFileWriter openSAMFileWriterToFile(SAMFileWriterFactory sfw,File file,SAMFileHeader header,boolean presorted);
public abstract boolean isBinary();

}
