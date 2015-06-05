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
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.knime.AbstractKnimeApplication;
import com.github.lindenb.jvarkit.util.illumina.FastQName;

public class IlluminaDirectory
	extends AbstractKnimeApplication
	{
	private int ID_GENERATOR=0;
	private boolean JSON=false;
    private MessageDigest md5;

    
    public void setJSON(boolean jSON) {
		JSON = jSON;
		}
    
    public boolean isJSON() {
		return JSON;
		}
    
    private  String md5(String in)
    	{
    	if(md5==null)
	    	{
	    	  try {
	              md5 = MessageDigest.getInstance("MD5");
	          } catch (NoSuchAlgorithmException e) {
	              throw new RuntimeException("MD5 algorithm not found", e);
	          }
	    	}
    	 md5.reset();
         md5.update(in.getBytes());
         String s = new BigInteger(1, md5.digest()).toString(16);
         if (s.length() != 32) {
             final String zeros = "00000000000000000000000000000000";
             s = zeros.substring(0, 32 - s.length()) + s;
         }
         return s;
    	}
    
    private class Folder
    	{
    	SortedMap<String, Sample> sampleMap=new TreeMap<String, IlluminaDirectory.Sample>();
    	List<Pair> undtermined=new ArrayList<Pair>();
    	void scan(File f)
    		{
    		if(f==null) return;
    		if(!f.canRead()) return;
    		IlluminaDirectory.this.info("Scanning "+f);
    		
			FastQName fq=FastQName.parse(f);
			if(!fq.isValid())
				{
				warning("invalid name:"+fq);
				return;
				}
			if(fq.isUndetermined())
				{
				for(int i=0;i< undtermined.size();++i)
	    			{
	    			Pair p=undtermined.get(i);
	    			if(p.complement(fq)) return;
	    			}
				undtermined.add(new Pair(fq));
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
				}
				
    			
    		}
    	
    	void json(PrintStream out)
    		{
    		boolean first=true;
    		out.print("{\"samples\":[");
    		for(Sample S:this.sampleMap.values())
				{
    			if(!first) out.print(',');
    			first=false;
				S.json(out);
				}
    		out.print("],\"undetermined\":[");
    		first=true;
    		for(Pair p:undtermined)
				{
    			if(!first) out.print(',');
    			first=false;
				p.json(out);
				}
    		out.print("]}");
    		}
    	
    	void write(XMLStreamWriter w) throws XMLStreamException
    		{
    		w.writeStartElement("project");
    		w.writeAttribute("name", "Project1");
    		w.writeAttribute("center", "CENTER");
    		w.writeAttribute("haloplex", "false");
    		w.writeAttribute("wgs", "false");

    		for(Sample S:this.sampleMap.values())
    			{
    			S.write(w);
    			}
    		w.writeStartElement("undetermined");
    		for(Pair p:undtermined)
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
    	
    	boolean complement(FastQName other)
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
    	
    	void json(PrintStream out)
    		{
    		if(forward!=null && reverse!=null)
				{
	    		
	    		out.print("{");
	    		out.print("\"id\":\"p"+this.id+"\",");
	    		out.print("\"md5pair\":\""+md5(forward.getFile().getPath()+reverse.getFile().getPath())+"\",");
	    		out.print("\"lane\":"+forward.getLane()+",");
    			if(forward.getSeqIndex()!=null)
    				{
    				out.print("\"index\":\""+forward.getSeqIndex()+"\",");
    				}
    			else
    				{
    				out.print("\"index\":null,");
    				}
	    		out.print("\"split\":"+forward.getSplit()+",");
	    		out.print("\"forward\":{");
	    		out.print("\"md5filename\":\""+md5(forward.getFile().getPath())+"\",");
	    		out.print("\"path\":\""+forward.getFile().getPath()+"\",");
	    		out.print("\"side\":"+forward.getSide().ordinal()+",");
	    		out.print("\"file-size\":"+forward.getFile().length());
	    		out.print("},\"reverse\":{");
	    		out.print("\"md5filename\":\""+md5(reverse.getFile().getPath())+"\",");
	    		out.print("\"path\":\""+reverse.getFile().getPath()+"\",");
	    		out.print("\"side\":"+reverse.getSide().ordinal()+",");
	    		out.print("\"file-size\":"+reverse.getFile().length());
	    		out.print("}}");
				}
    		else
    			{
    			FastQName F=(forward==null?reverse:forward);
    			out.print("{");
    			out.print("\"id\":\"p"+this.id+"\",");
    			out.print("\"md5filename\":\""+md5(F.getFile().getPath())+"\",");
	    		out.print("\"lane\":"+F.getLane()+",");
    			if(forward.getSeqIndex()!=null)
    				{
    				out.print("\"index\":\""+F.getSeqIndex()+"\",");
    				}
    			else
    				{
    				out.print("\"index\":null,");
    				}
	    		out.print("\"split\":"+F.getSplit()+",");
	    		out.print("\"path\":\""+F.getFile().getPath()+"\",");
	    		out.print("\"side\":"+F.getSide()+",");
    			out.print("}");
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
    	}
    
    private class Sample
		{
		String name;
    	List<Pair> pairs=new ArrayList<Pair>();
    	
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
    	
    	void json(PrintStream out)
    		{
    		out.print("{\"sample\":\""+ this.name +"\",\"files\":[");
    		for(int i=0;i< pairs.size();++i)
    			{
    			if(i>0) out.print(",");
    			pairs.get(i).json(out);
    			}
    		out.print("]}");
    		}
    	
    	
		}
    
    
    
    @Override
	public String getProgramDescription() {
		return "Scan folders and generate a structured summary of the files in JSON or XML";
		}
    @Override
    protected String getOnlineDocUrl() {
    	return DEFAULT_WIKI_PREFIX+"Illuminadir";
    	}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -J produces a JSON output");
		super.printOptions(out);
		}
	

	@Override
	public int doWork(String[] args)
			{
			com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
			int c;
			while((c=opt.getopt(args, super.getGetOptDefault()+"Jo:"))!=-1)
				{
				switch(c)
					{
					case 'o': setOutputFile(opt.getOptArg()); break;
					case 'J': setJSON(true);break;
					default:
						{
						switch(super.handleOtherOptions(c, opt, args))
							{
							case EXIT_FAILURE:return -1;
							case EXIT_SUCCESS:return 0;
							case OK:break;
							}
						}
					}
				}
			return mainWork(opt.getOptInd(), args);
			}
	
	
		@Override
		public int executeKnime(List<String> args)
			{
			if(!args.isEmpty())
				{
				error("Expected to read filenames on stdin.");
				return -1;
				}
			
			try
				{
				List<Folder> folders=new ArrayList<Folder>();
				Folder folder=new Folder();
				folders.add(folder);
				String line;
				BufferedReader in=new BufferedReader(new InputStreamReader(System.in));
				while((line=in.readLine())!=null)
					{
					if(line.isEmpty() || line.startsWith("#")) continue;
					if(!line.endsWith(".fastq.gz"))
						{
						warning("ignoring "+line);
						continue;
						}
					File f=new File(line);
					if(!f.exists())
						{
						error("Doesn't exist:"+f);
						return -1;
						}
					if(!f.isFile())
						{
						error("Not a file:"+f);
						return -1;
						}
					folder.scan(f);
					}
				in.close();
		    	
				PrintStream pw = getOutputFile()==null?System.out:new PrintStream(getOutputFile());

		    	if(this.isJSON())
		    		{
		    		pw.print("[");
		    		for(int i=0;i< folders.size();++i)
		    			{
		    			if(i>0) pw.print(",");
		    			folders.get(i).json(pw);
		    			}
		    		
		    		pw.println("]");
		    		pw.flush();
		    		int ret = pw.checkError()?-1:0;
		    		pw.close();
		    		return ret;
		    		}
		    	else
		    		{
	    			XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
	    			XMLStreamWriter w= xmlfactory.createXMLStreamWriter(pw);
	    			w.writeStartDocument("UTF-8","1.0");
	    			w.writeStartElement("model");
	    			w.writeComment(this.getProgramCommandLine());
	    			for(Folder f:folders) f.write(w);
	    			w.writeEndElement();
	    			w.writeEndDocument();
	    			w.flush();
	    			w.close();
		    		}
				pw.close();
				}
			catch(Exception err)
				{
				error(err);
				return -1;
				}
			finally
				{
				
				}
			return 0;
			}
   

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new IlluminaDirectory().instanceMainWithExit(args);
		}

}
