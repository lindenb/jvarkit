/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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
* 2017 creation

*/
package com.github.lindenb.jvarkit.util.bio.fasta;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;

public abstract class ContigNameConverter implements  Function<String, String> {
public static final String OPT_ON_NT_FOUND_DESC="Contig converter. I will do my best to convert the contig names (e.g 'chr1' -> '1'): But what should I do when comparing two dictionaries with different notations";	
public static enum OnNotFound{RAISE_EXCEPTION,SKIP,RETURN_ORIGINAL};
private OnNotFound onNotFound=OnNotFound.RAISE_EXCEPTION;
protected abstract String find(final String contig);

public String getName() {
	return "Untitled";
}

public ContigNameConverter  setOnNotFound(final OnNotFound onNotFound) {
	this.onNotFound = onNotFound;
	return this;
	}

public OnNotFound getOnNotFound() {
	return onNotFound;
}

@Override
public String apply(final String contig) {
	final String dest=find(contig);
	if(dest==null)
		{
		switch(onNotFound)
			{
			case RAISE_EXCEPTION:
				throw new JvarkitException.ContigNotFound("cannot convert contig \""+contig+"\"");
			case SKIP: /* lave plus blanc */
				return null;
			case RETURN_ORIGINAL:
				return contig;
			}
		}
	return dest;
	}		

public SAMSequenceRecord convertSAMSequenceRecord(final SAMSequenceRecord ssr) {
	final String newName=apply(ssr.getSequenceName());
	if(newName==null)
		{
		return null;
		}
	final SAMSequenceRecord ssr2=new SAMSequenceRecord(newName, ssr.getSequenceLength());
		
	return ssr2;
	}


public SAMSequenceDictionary convertDictionary(final SAMSequenceDictionary dict) {
	return new SAMSequenceDictionary(dict.getSequences().
		stream().
		map(S->convertSAMSequenceRecord(S)).
		filter(S->S!=null).
		collect(Collectors.toList())
		);
	}

public static  ContigNameConverter getIdentity() {
	return new ContigNameConverter()
			{
			@Override
			protected String find(final String contig) {
				return contig;
				}
			@Override
			public String getName() {
				return "Identity";
			}
			};
	}

private static class MapBasedContigNameConverter extends ContigNameConverter
	{
	protected String name="Untitled";
	protected final Map<String,String> map=new HashMap<>();
	@Override
	protected String find(final String contig) {
		return this.map.get(contig);
		}
	@Override
	public String getName() {
		return name;
		}
	
	}

public static ContigNameConverter fromFile(final File mappingFile)
	{
	IOUtil.assertFileIsReadable(mappingFile);
	final MapBasedContigNameConverter mapper=new MapBasedContigNameConverter();
	mapper.name=mappingFile.getName();
	BufferedReader in=null;
	try {
		in=IOUtils.openFileForBufferedReading(mappingFile);
		String line;
		while((line=in.readLine())!=null)
			{
			if(line.isEmpty() || line.startsWith("#")) continue;
			final String tokens[]=line.split("[\t]");
			if(tokens.length!=2
					|| tokens[0].trim().isEmpty()
					|| tokens[1].trim().isEmpty()
					) {
				in.close();in=null;
				throw new IOException("Bad mapping line: \""+line+"\"");
				}
			tokens[0]=tokens[0].trim();
			tokens[1]=tokens[1].trim();
			if(mapper.map.containsKey(tokens[0]))
				{
				in.close();
				throw new IOException("Mapping defined twice for: \""+tokens[0]+"\"");
				}
			mapper.map.put(tokens[0], tokens[1]);
			}
		return mapper;
		}
	catch(final IOException err) {
		throw new RuntimeIOException(err);
		}
	finally {
		CloserUtil.close(in);
	}
	
}

private static class GRCh37Ucsc2Ensembl extends MapBasedContigNameConverter
	{
	GRCh37Ucsc2Ensembl() {
		for(int i=1;i<=22;++i) map.put("chr"+i, String.valueOf(i));
		map.put("chrX","X");map.put("chrY","Y");map.put("chrM","MT");
		}
	@Override
	public String getName() {return "GRCh37UcscToEnsembl";}
	}

/** creates a ContigNameConverter from two dictionaries, just using the common chromosome names */
public static ContigNameConverter fromDictionaries(
		final SAMSequenceDictionary dictIn,
		final SAMSequenceDictionary dictOut
		)
	{
	return new TwoDictionaries(dictIn,dictOut);
	}

/** creates a ContigNameConverter from one dictionary, try to guess best conversion*/
public static ContigNameConverter fromOneDictionary(
		final SAMSequenceDictionary dict
		)
	{
	return new OneDictionary(dict);
	}
/** creates a ContigNameConverter from a set of contig names, try to guess best conversion*/
public static ContigNameConverter fromContigSet(
		final Set<String> contigSet
		)
	{
	return new OneDictionary(contigSet);
	}



private static class TwoDictionaries extends MapBasedContigNameConverter
	{
	TwoDictionaries(final SAMSequenceDictionary dictIn,final SAMSequenceDictionary dictOut) {
		if(dictIn==null || dictIn.isEmpty())
			{
			throw new IllegalArgumentException("dictIn is null or empty");
			}
		if(dictOut==null || dictOut.isEmpty())
			{
			throw new IllegalArgumentException("dictOut is null or empty");
			}
		
		for(final SAMSequenceRecord ssr : dictIn.getSequences())
			{
			if(dictOut.getSequence(ssr.getSequenceName())!=null)
				{
				super.map.put(ssr.getSequenceName(), ssr.getSequenceName());
				continue;
				}
			
			if(ssr.getSequenceName().startsWith("chr") )
				{
				final String ctg = ssr.getSequenceName().substring(3);
				if(dictOut.getSequence(ctg)!=null)
					{
					super.map.put(ssr.getSequenceName(), ctg);
					continue;
					}
				if(ssr.getSequenceName().equals("chrMT") && dictOut.getSequence("M")!=null)
					{
					super.map.put("chrMT","M");
					continue;
					}
				if(ssr.getSequenceName().equals("chrM") && dictOut.getSequence("MT")!=null)
					{
					super.map.put("chrM","MT");
					continue;
					}
				continue;
				}
			if(!ssr.getSequenceName().startsWith("chr") )
				{
				final String ctg = "chr"+ssr.getSequenceName();
				if(dictOut.getSequence(ctg)!=null)
					{
					super.map.put(ssr.getSequenceName(), ctg);
					continue;
					}
				if(ssr.getSequenceName().equals("MT") && dictOut.getSequence("chrM")!=null)
					{
					super.map.put("MT","chrM");
					continue;
					}
				if(ssr.getSequenceName().equals("M") && dictOut.getSequence("chrMT")!=null)
					{
					super.map.put("M","chrMT");
					continue;
					}
				continue;
				}
			
			final String md5in = ssr.getMd5();
			if(!StringUtil.isBlank(md5in)) {
				final List<SAMSequenceRecord> ssrOUts = dictOut.getSequences().
						stream().
						filter(SSR->!StringUtil.isBlank(SSR.getMd5()) && md5in.equals(SSR.getMd5())).
						collect(Collectors.toList());
				if(ssrOUts.size()==1)
					{
					super.map.put(ssr.getSequenceName(), ssrOUts.get(0).getSequenceName());
					continue;
					}
				}
			}
		}
	@Override
	public String getName() {return "TwoDictionaries";}
	}

private static class OneDictionary extends ContigNameConverter
	{
	private final SAMSequenceDictionary dict;
	OneDictionary( final SAMSequenceDictionary dict)
		{
		this.dict = dict;
		}
	
	OneDictionary( final Set<String> contigs)
		{
		this.dict = new SAMSequenceDictionary(
			contigs.stream().
				map(C->new SAMSequenceRecord(C, 9999)).
				collect(Collectors.toList())
			);
		}
	
	@Override
	protected String find(final String contig) {
		if(this.dict.getSequenceIndex(contig)!=-1) return contig;
		
	
		if(contig.equals("chrM"))
			{
			SAMSequenceRecord ssr = this.dict.getSequence("M");
			if(ssr==null) ssr = this.dict.getSequence("MT");
			if(ssr!=null) return ssr.getSequenceName();
			}
		else if(contig.equals("M") || contig.equals("MT"))
			{
			SAMSequenceRecord ssr = this.dict.getSequence("chrM");
			if(ssr==null) ssr = this.dict.getSequence("chrMT");
			if(ssr!=null) return ssr.getSequenceName();
			}
		
		final String c2;
		if(contig.startsWith("chr"))
			{
			c2 = contig.substring(3);
			}
		else
			{
			c2 = "chr"+contig;
			}
		final SAMSequenceRecord ssr = this.dict.getSequence(c2);
		if(ssr!=null) return ssr.getSequenceName();
		
		return null;
		}
	@Override
	public String getName() {
		return "OneDictionary";
		}
	}

}
