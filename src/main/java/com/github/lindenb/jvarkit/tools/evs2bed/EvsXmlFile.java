/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.evs2bed;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.TribbleIndexedFeatureReader;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.readers.LineIterator;
public class EvsXmlFile implements Closeable{
private Index tribbleIndex=null;
private TribbleIndexedFeatureReader<SnpDataFeature, LineIterator> tribble;

/** 
 * tribble access to XML file generated with evsDumpXml
 * @param xmlFile
 * @throws IOException
 */
public EvsXmlFile(File xmlFile) throws IOException
	{	
	File idxFile = Tribble.indexFile(xmlFile);
	SnpDataCodec tribbleCodec=new SnpDataCodec();
	if(idxFile.exists() && idxFile.lastModified()>=xmlFile.lastModified())
	    {
		tribbleIndex=IndexFactory.loadIndex(idxFile.getPath());
	    }
	else
	    {
	 	tribbleIndex=IndexFactory.createLinearIndex(xmlFile, tribbleCodec);
	    }
	
	this.tribble =
			new TribbleIndexedFeatureReader<>(
					xmlFile.getPath(),
					tribbleCodec,
					tribbleIndex
					);

	}

@Override
protected void finalize() throws Throwable {
	this.close();
	super.finalize();
	}

@Override
public void close() throws IOException {
	if(this.tribble!=null) this.tribble.close();
	this.tribble=null;
	this.tribbleIndex=null;
	}

public CloseableIterator<SnpDataFeature> iterator()
	throws IOException
	{
	if(this.tribble==null) throw new java.io.IOException("Tribble file was closed");
	return this.tribble.iterator();
	}
public CloseableIterator<SnpDataFeature> query(String chrom,int start,int end)
	throws IOException
	{
	if(this.tribble==null) throw new java.io.IOException("Tribble file was closed");
	return this.tribble.query(chrom,start,end);
	}
}
