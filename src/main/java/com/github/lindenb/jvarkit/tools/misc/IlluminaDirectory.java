/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.misc;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.illumina.FastQName;
import com.google.gson.stream.JsonWriter;

import htsjdk.samtools.util.CloserUtil;

public class IlluminaDirectory
	extends AbstractIlluminaDirectory
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(IlluminaDirectory.class);

	private int ID_GENERATOR=0;

    
    private class Folder
    	{
    	String projectName="Project1";
    	final SortedMap<String, Sample> sampleMap=new TreeMap<String, IlluminaDirectory.Sample>();
    	final List<Pair> undetermined=new ArrayList<Pair>();
    	void scan(final File f)
    		{
    		if(f==null) return;
    		if(!f.canRead()) return;
    		LOG.info("Scanning "+f);
    		
			FastQName fq=FastQName.parse(f);
			if(!fq.isValid())
				{
				LOG.warn("invalid name:"+fq);
				return;
				}
			if(fq.isUndetermined())
				{
				for(int i=0;i< undetermined.size();++i)
	    			{
	    			Pair p=undetermined.get(i);
	    			if(p.complement(fq)) return;
	    			}
				undetermined.add(new Pair(fq));
				}
			else
				{
				Sample sample=this.sampleMap.get(fq.getSample());
				if(sample==null)
					{
					sample=new Sample();
					sample.name=fq.getSample();
					this.sampleMap.put(sample.name,sample);
					}
				sample.add(fq);
				
				
				final File sampleDir = f.getParentFile();
				if(sampleDir!=null && sampleDir.isDirectory() && sampleDir.getName().startsWith("Sample_"))
					{
					final File projDir = sampleDir.getParentFile();
					if(projDir!=null && projDir.isDirectory() && projDir.getName().startsWith("Project_"))
						{
						this.projectName = projDir.getName().substring(8).replace(' ', '_');
						}
					}
				}
				
    			
    		}
    	
    	void json(final JsonWriter out)  throws IOException
    		{
    		out.beginObject();
    		out.name("samples");
    		out.beginArray();

    		for(final Sample S:this.sampleMap.values())
				{
				S.json(out);
				}
    		out.endArray();
    		out.name("undetermined");
    		out.beginArray();
    		for(final Pair p:undetermined)
				{
				p.json(out);
				}
    		out.endArray();
    		out.endObject();
    		}
    	
    	void write(final XMLStreamWriter w) throws XMLStreamException
    		{
    		w.writeStartElement("project");
    		w.writeAttribute("name",this.projectName);
    		w.writeAttribute("center", "CENTER");
    		w.writeAttribute("haloplex", "false");
    		w.writeAttribute("wgs", "false");

    		for(final Sample S:this.sampleMap.values())
    			{
    			S.write(w);
    			}
    		w.writeStartElement("undetermined");
    		for(final Pair p:this.undetermined)
    			{
    			p.write(w);
    			}
    		w.writeEndElement();
    		w.writeEndElement();
    		}
    	
    	}
    
    /** 
     * A pair of fastq , Forward, reverse
     */
    private class Pair
    	{
    	int id;
    	FastQName forward;
    	FastQName reverse;
    	
    	Pair(FastQName fq)
    		{
    		id=++ID_GENERATOR;
    		switch(fq.getSide())
    			{
    			case Forward:forward=fq; break;
    			case Reverse:reverse=fq; break;
    			default:throw new RuntimeException("bad side "+fq);
    			}
    		}
    	
    	boolean complement(final FastQName other)
    		{
    		if(forward!=null && reverse!=null) return false;
    		if(forward!=null && forward.isComplementOf(other))
    			{
    			reverse=other;
    			return true;
    			}
    		else if(reverse!=null && reverse.isComplementOf(other))
    			{
    			forward=other;
    			return true;
    			}
    		return false;
    		}
    	
    	void json(final JsonWriter out) throws IOException
    		{
    		if(forward!=null && reverse!=null)
				{
	    		out.beginObject();
	    		out.name("id");out.value("p"+this.id);
	    		out.name("md5pair");out.value(md5(forward.getFile().getPath()+reverse.getFile().getPath()));
	    		out.name("lane");out.value(""+forward.getLane());
	    		out.name("index");
	    		if(forward.getSeqIndex()!=null)
    				{
    				out.value(forward.getSeqIndex());
    				}
    			else
    				{
    				out.nullValue();
    				}
	    		out.name("split");out.value(""+forward.getSplit());
	    		
	    		
	    		out.name("forward");
	    		
	    		out.beginObject();
	    		out.name("md5filename");out.value(md5(forward.getFile().getPath()));
	    		out.name("path");out.value(forward.getFile().getPath());
	    		out.name("side");out.value(forward.getSide().ordinal());
	    		out.name("file-size");out.value(forward.getFile().length());	    		
	    		out.endObject();
	    		
	    		out.name("reverse");
	    		
	    		out.beginObject();
	    		out.name("md5filename");out.value(md5(reverse.getFile().getPath()));
	    		out.name("path");out.value(reverse.getFile().getPath());
	    		out.name("side");out.value(reverse.getSide().ordinal());
	    		out.name("file-size");out.value(reverse.getFile().length());	    		
	    		out.endObject();

	    		out.endObject();
				}
    		else
    			{
    			final FastQName F=(forward==null?reverse:forward);
    			out.beginObject();
	    		out.name("id");out.value("p"+this.id);
	    		out.name("md5filename");out.value(md5(F.getFile().getPath()));
	    		out.name("lane");out.value(""+F.getLane());
	    		out.name("index");
	    		if(F.getSeqIndex()!=null)
    				{
    				out.value(F.getSeqIndex());
    				}
    			else
    				{
    				out.nullValue();
    				}
	    		out.name("split");out.value(""+F.getSplit());
	    		out.name("path");out.value(F.getFile().getPath());
	    		out.name("side");out.value(F.getSide().ordinal());

	    		
    			out.endObject();
    			}
			}
    	
    	void write(XMLStreamWriter w,String tagName,FastQName fastqFile) throws XMLStreamException
    		{
			w.writeStartElement(tagName);
			w.writeAttribute("md5filename",md5(fastqFile.getFile().getPath()));
			w.writeAttribute("file-size",String.valueOf( fastqFile.getFile().length()));
			w.writeCharacters(fastqFile.getFile().getPath());
			w.writeEndElement();
    		}
    	
    	void write(XMLStreamWriter w) throws XMLStreamException
    		{
			w.writeStartElement("fastq");
			w.writeAttribute("id","p"+this.id);
			w.writeAttribute("md5",md5(forward.getFile().getPath()+reverse.getFile().getPath()));
			w.writeAttribute("lane", String.valueOf(forward.getLane()));
			if(forward.getSeqIndex()!=null) w.writeAttribute("index", String.valueOf(forward.getSeqIndex()));
			w.writeAttribute("split", String.valueOf(forward.getSplit()));
			w.writeAttribute("group-id", getGroupId());
			
			if(forward!=null && reverse!=null)
    			{
    			write(w,"for",forward);
    			write(w,"rev",reverse);
    			}
			else
				{
				write(w,"single",forward==null?reverse:forward);
				}
			
		
			w.writeEndElement();
    		}
    	private String getGroupId()
    		{
    		return IlluminaDirectory.this.getGroupId(this.forward);
    		}
    	}
    
    private Map<String,String> groupIdMap=new HashMap<>(); 
    private String getGroupId(final FastQName fastq)
    	{
    	final String s= md5( fastq.getSample()+" "+fastq.getLane()+" "+fastq.getSeqIndex());
    	String gid = this.groupIdMap.get(s);
    	if(gid==null)
    		{
    		gid = fastq.getSample()+"."+(this.groupIdMap.size()+1);
    		this.groupIdMap.put(s, gid);
    		}
    	return gid;
    	}
    
    
    private class Sample
		{
		String name;
    	final List<Pair> pairs=new ArrayList<Pair>();
    	
    	private void add(FastQName fq)
    		{
    		for(int i=0;i< pairs.size();++i)
    			{
    			Pair p=pairs.get(i);
    			if(p.complement(fq)) return;
    			}
    		pairs.add(new Pair(fq));
    		}
    	
    	void write(XMLStreamWriter w) throws XMLStreamException
			{
			w.writeStartElement("sample");
			w.writeAttribute("name",this.name);
			w.writeAttribute("father","undefined");
			w.writeAttribute("mother","undefined");
			w.writeAttribute("sex","undefined");
			
			for(Pair p:this.pairs)
				{
				p.write(w);
				}
				
			w.writeEndElement();
			}
    	
    	void json(final JsonWriter out)  throws IOException
    		{
    		out.beginObject();
    		out.name("sample");
    		out.value(this.name);
    		out.name("files");
    		out.beginArray();
    		
    		for(final Pair p: this.pairs)
    			{
    			p.json(out);
    			}
    		
    		out.endArray();
    		out.endObject();
    		}
    	
    	
		}
    
	    @Override
	    protected Collection<Throwable> call(final String inputName) throws Exception {
	    	BufferedReader in=null;;
			try
				{
				if( inputName == null) 
					{
					in = IOUtils.openStreamForBufferedReader(stdin());
					}
				else
					{
					in = IOUtils.openURIForBufferedReading(inputName);
					}
				
				final Folder folder=new Folder();
				String line;
				while((line=in.readLine())!=null)
					{
					if(line.isEmpty() || line.startsWith("#")) continue;
					if(!line.endsWith(".fastq.gz"))
						{
						LOG.warn("ignoring "+line+" because it doesn't end with *.fastq.gz");
						continue;
						}
					final File f=new File(line);
					if(!f.exists())
						{
						return wrapException("Doesn't exist:"+f);
						}
					if(!f.isFile())
						{
						return wrapException("Not a file:"+f);
						}
					folder.scan(f);
					}
				in.close();
		    	
				final PrintWriter pw = this.openFileOrStdoutAsPrintWriter();

		    	if(super.isJSON())
		    		{
		    		final JsonWriter js=new JsonWriter(pw);
		    		folder.json(js);
		    		CloserUtil.close(js);
		    		}
		    	else
		    		{
		    		final XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
		    		final XMLStreamWriter w= xmlfactory.createXMLStreamWriter(pw);
	    			w.writeStartDocument("UTF-8","1.0");
	    			folder.write(w);
	    			w.writeEndDocument();
	    			w.flush();
	    			CloserUtil.close(w);
		    		}
		    	pw.flush();
		    	CloserUtil.close(pw);
		    	return RETURN_OK;
				}
			catch(final Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(in);
				}
			}
   

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new IlluminaDirectory().instanceMainWithExit(args);
		}

}
