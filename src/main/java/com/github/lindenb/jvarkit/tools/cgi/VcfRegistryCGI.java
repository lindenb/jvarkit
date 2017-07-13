package com.github.lindenb.jvarkit.tools.cgi;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
/**
BEGIN_DOC

## Motivation

CGI/Web tool printing all the variants at a given position for a collection of tabix-ed VCF files.

## output example:
<table><thead><tr><th>CHROM</th><th>POS</th><th>ID</th><th>REF</th><th>QUAL</th><th>Sample</th><th>Alleles</th><th>DP</th><th>GQ</th><th>File</th></tr></thead><tbody><tr><td>22</td><td>11905981</td><td><a href="http://www.ncbi.nlm.nih.gov/snp/86876">rs86876</a></td><td>A</td><td>3133</td><td>SAMPLE_273</td><td><span style="color:blue;">A/G</span></td><td>211</td><td>99</td><td>/path/to/SAMPLE_273_variations.gatk.annotations.vcf.gz</td></tr><tr><td>22</td><td>11905981</td><td><a href="http://www.ncbi.nlm.nih.gov/snp/86876">rs86876</a></td><td>A</td><td>143</td><td>SAMPLE_273</td><td><span style="color:red;">G/G</span></td><td>221</td><td>99</td><td>/path/to/SAMPLE_273_variations.samtools.annotations.vcf.gz</td></tr><tr><td>22</td><td>11905981</td><td><a href="http://www.ncbi.nlm.nih.gov/snp/86876">rs86876</a></td><td>A</td><td>2455</td><td>SAMPLE_324</td><td><span style="color:blue;">A/G</span></td><td>165</td><td>99</td><td>/path/to/SAMPLE_324_variations.gatk.annotations.vcf.gz</td></tr><tr><td>22</td><td>11905981</td><td><a href="http://www.ncbi.nlm.nih.gov/snp/86876">rs86876</a></td><td>A</td><td>118</td><td>SAMPLE_324</td><td><span style="color:red;">G/G</span></td><td>167</td><td>99</td><td>/path/to/SAMPLE_324_variations.samtools.annotations.vcf.gz</td></tr><tr><td>22</td><td>11905981</td><td><a href="http://www.ncbi.nlm.nih.gov/snp/86876">rs86876</a></td><td>A</td><td>6398</td><td>SAMPLE_016</td><td><span style="color:blue;">A/G</span></td><td>462</td><td>99</td><td>/path/to/SAMPLE_016_variations.gatk.annotations.vcf.gz</td></tr><tr><td>22</td><td>11905981</td><td><a href="http://www.ncbi.nlm.nih.gov/snp/86876">rs86876</a></td><td>A</td><td>133</td><td>SAMPLE_016</td><td><span style="color:blue;">A/G</span></td><td>482</td><td>99</td><td>/path/to/SAMPLE_016_variations.samtools.annotations.vcf.gz</td></tr> <tr><td>(...)</td><td>(...)</td><td><a (...)</td><td>(...)</td><td>(...)</td><td>(...)</td><td>(...)</td><td>(...)</td><td>(...)</td><td>(...)</td></tr></table>

##Compilation/Installation

### Create a Preference file.
create/update a java preferences file with the following key/example:

```xml
<?xml version="1.0"?>
<!DOCTYPE preferences SYSTEM "http://java.sun.com/dtd/preferences.dtd">
<preferences EXTERNAL_XML_VERSION="1.0">
  <root type="user">
    <map>
      <entry key="vcf.registry.group.file" value="/var/www/cgi-bin/vcfregistry.tsv"/>
    </map>
  </root>
</preferences>
```
here /var/www/cgi-bin/vcfregistry.tsv  is a tab delimited file containing the path to some collections of VCFs (path/Title):

```
$ cat vcfregistry.tsv
/var/www/cgi-bin/vcfregistry1.tsv       Collection1
/var/www/cgi-bin/vcfregistry2.tsv       Collection2
/var/www/cgi-bin/vcfregistry3.tsv       Collection3
```

A collection file contains a list of VCF indexed with tabix:

```
$ head /var/www/cgi-bin/vcfregistry1.tsv 
/commun/data/projects/TEST/align/Samples/SAMPLE_73/VCF/SAMPLE_73_variations.gatk.annotations.vcf.gz
/commun/data/projects/TEST/align/Samples/SAMPLE_73/VCF/SAMPLE_73_variations.samtools.annotations.vcf.gz
/commun/data/projects/TEST/align/Samples/SAMPLE_72/VCF/SAMPLE_72_variations.gatk.annotations.vcf.gz
/commun/data/projects/TEST/align/Samples/SAMPLE_72/VCF/SAMPLE_72_variations.samtools.annotations.vcf.gz
/commun/data/projects/TEST/align/Samples/SAMPLE_71/VCF/SAMPLE_71_variations.gatk.annotations.vcf.gz
/commun/data/projects/TEST/align/Samples/SAMPLE_71/VCF/SAMPLE_71_variations.samtools.annotations.vcf.gz
/commun/data/projects/TEST/align/Samples/SAMPLE_24/VCF/SAMPLE_24_variations.gatk.annotations.vcf.gz
/commun/data/projects/TEST/align/Samples/SAMPLE_24/VCF/SAMPLE_24_variations.samtools.annotations.vcf.gz
/commun/data/projects/TEST/align/Samples/SAMPLE_16/VCF/SAMPLE_16_variations.gatk.annotations.vcf.gz
/commun/data/projects/TEST/align/Samples/SAMPLE_16/VCF/SAMPLE_16_variations.samtools.annotations.vcf.gz
```

### Compile the program

See also [[Compilation]].

Important : the JAR/library files of picard should be visible from the CGI-bin. The path defined in the build.properties should be absolute.

```bash
$ make vcfregistry.cgi
```
It creates a jar file (vcfregistry.jar) and an executable shell script (vcfregistry.cgi) in the 'dist' folder. Both should be moved in the cgi-bin script of your server:
```bash
$ sudo mv dist/vcfregistry* /var/www/cgi-bin/
```

If needed, edit the script tviewweb.cgi and change the JVM property `Dprefs.file.xml=` to the correct place of your xml preference file.
```
java (...) -Dprefs.file.xml=/var/www/cgi-bin/prefs.xml (...)
````

END_DOC
 */
@Program(name="vcfregistry.cg",description="CGI/Web tool printing all the variants at a given position for a collection of tabix-ed VCF files.")
public class VcfRegistryCGI extends AbstractCGI {
	private static final Logger LOG=Logger.build(VcfRegistryCGI.class).make();

	 private static final String GROUPID_PARAM="g";
	    private static final String RGN_PARAM="r";
	    private Throwable lastException=null;
	    private static class Position
	    	{
	    	String chrom;
	    	int pos;
	    	}
	    
	    private static class FileWithComment
	    	{
	    	File file;
		    String _desc;
		    FileWithComment(File file)
		     	{
		    	this.file=file;
		    	this._desc="";
		     	}
		    public String getDesc()
		    	{
		    	return _desc==null || _desc.isEmpty()?file.getName():_desc;
		    	}
	    	}
	    
	    /** list of VCF path stored in a file */
	    private class GroupFile
	    	extends FileWithComment
	        {
	    	int index=0;
	    	GroupFile(File file)
	    		{
	    		super(file);
	    		}
	        }
	    
	    private class VcfFile
    	extends FileWithComment
	        {
	    	VcfFile(File file)
	    		{
	    		super(file);
	    		}
	        }
	   
	    private VcfRegistryCGI()
	        {
	        }
	   
	    
	    private Position parsePosition()
	    	{
	    	return parsePosition(getString(RGN_PARAM));
	    	}
	    
	    private Position parsePosition(String s)
	    	{
	    	if(s==null) return null;
	    	int colon=s.indexOf(':');
	    	if(colon<1) return null;
	    	int pos;
	    	try
	    		{
	    		pos=Integer.parseInt(s.substring(colon+1).trim());
	    		}
	    	catch(Exception err)
	    		{
	    		return null;
	    		}
	    	if(pos<1) return null;
	    	Position position=new Position();
	    	position.chrom=s.substring(0, colon).trim();
	    	if(position.chrom.trim().isEmpty()) return null;
	    	position.pos=pos;
	    	return position;
	    	}	
	    
	    /* get the path to the main config file */
	    private File getGroupFile()
	        {
	        String groupFileStr=null;
	        try
	            {
	            groupFileStr=getPreferences().get("vcf.registry.group.file", null);
	            if(groupFileStr==null) return null;
	            File f= new File(groupFileStr);
	            if(f.isFile() && f.canRead()) return f;
	            return null;
	            }
	        catch(IOException err)
	            {
	        	lastException=err;
	            return null;
	            }
	        }
	    
	    /* read the main config file and return all the group files */
	    private List<GroupFile> getGroupFiles()
	        {
	        File g=getGroupFile();
	        if(g==null) return Collections.emptyList();
	        BufferedReader in=null;
	        List<GroupFile> L=new ArrayList<GroupFile>();
	        try
	            {
	            in=IOUtils.openFileForBufferedReading(g);
	            int tab;
	            String line;
	            while((line=in.readLine())!=null)
	                {
	                if(line.startsWith("#") || line.isEmpty() || (tab=line.indexOf('\t'))<1)
	                    continue;
	                File f=new File(line.substring(0,tab).trim());
	                GroupFile gf=new GroupFile(f);
	                if(!(gf.file.exists() && gf.file.isFile() && gf.file.canRead()))
	                		{
	                		in.close();
	                		in=null;
	                		throw new IOException("Bad file in "+line+" / "+g);
	                		}
	                gf._desc=line.substring(tab+1).trim();
	                gf.index=L.size();
	                L.add(gf);
	                }
	            return L;
	            }
	        catch(Exception err)
	            {
	        	this.lastException=err;
	            return Collections.emptyList();
	            }
	        finally
	            {
	            CloserUtil.close(in);
	            }
	        }
	   
	    /** return all the VCF in a group file */
	    private List<VcfFile> getVcfFiles(GroupFile gf)
	        {
	        if(gf==null || gf.file==null) return Collections.emptyList();
	        BufferedReader in=null;
	        List<VcfFile> L=new ArrayList<VcfFile>();
	        try
	            {
	            in=IOUtils.openFileForBufferedReading(gf.file);
	            String line;
	            while((line=in.readLine())!=null)
	                {
	                if(line.startsWith("#") || line.isEmpty()) continue;
	                int tab=line.indexOf('\t');
	                if(tab==0) continue;
	                File f=new File(tab==-1?line:line.substring(0, tab));
	                if(!(f.exists() && f.isFile() && f.canRead()) )
	                	{
	                	in.close();in=null;
	                	throw new FileNotFoundException("Error for  "+f);
	                	}
	                File tbi=new File(line+".tbi");
	                if(!(tbi.exists() && tbi.isFile() && tbi.canRead()) )
	                	{
	                	continue;
	                	}
	                VcfFile vcf=new VcfFile(f);
	                if(tab>0) vcf._desc=line.substring(tab+1).trim();
	                L.add(vcf);
	                }
	            return L;
	            }
	        catch(IOException err)
	            {
	        	this.lastException=err;
	            return Collections.emptyList();
	            }
	        finally
	            {
	            CloserUtil.close(in);
	            }
	        }

	   
	   
	    private void welcomePane()
	        {
	        setMimeHeaderPrinted(true);
	        System.out.print("Content-type: text/html;charset=utf-8\n");
	        System.out.println();
	        System.out.flush();
	       
	        XMLStreamWriter w=null;
	        try
	            {
	            XMLOutputFactory xof=XMLOutputFactory.newFactory();
	            w=xof.createXMLStreamWriter(System.out,"UTF-8");
	            w.writeStartElement("html");
	            w.writeStartElement("body");
	           
	            w.writeStartElement("h1");
	            w.writeCharacters("VCFRegistry");
	            w.writeEndElement();
	           
	           
	            w.writeStartElement("ul");
	            for(GroupFile gf:getGroupFiles())
	                {
	                w.writeStartElement("li");
	                w.writeStartElement("a");
	                w.writeAttribute("title",gf.file.getPath());
	                w.writeAttribute("href","?"+GROUPID_PARAM+"="+gf.index);
	                w.writeCharacters(gf.getDesc());
	                w.writeEndElement();
	                w.writeEndElement();
	                }
	            w.writeEndElement();
	           
	            writeHTMLException(w, this.lastException);
	            
	            writeHtmlFooter(w);
	            w.writeEndElement();//body
	            w.writeEndElement();//html
	            }
	        catch(Exception err)
	            {
	        	
	            LOG.error(err);
	            }
	        finally
	            {
	            if(w!=null)
	                {
	                try {w.flush();} catch(XMLStreamException err){}
	                CloserUtil.close(w);
	                }
	            CloserUtil.close(System.out);
	            }
	        }
	   
	    private void printForm(
	    		XMLStreamWriter w,
	    		final GroupFile gf
	    		) throws XMLStreamException
	    	{
	    	 w.writeStartElement("form");
	    	
	    	 w.writeEmptyElement("input");
	    	 w.writeAttribute("type", "hidden");
	    	 w.writeAttribute("name",GROUPID_PARAM);
	    	 w.writeAttribute("value",String.valueOf(gf.index));
	    	 
	    	 w.writeStartElement("div");
	    	 
	    	 w.writeStartElement("label");
	    	 w.writeAttribute("for", "position");
	    	 w.writeCharacters("Position:");
	    	 w.writeEndElement();//label
	    	 
	    	 w.writeEmptyElement("input");
	    	 w.writeAttribute("id", "position");
	    	 w.writeAttribute("type", "text");
	    	 w.writeAttribute("name",RGN_PARAM);
	    	 w.writeAttribute("placeholder","chrom:position");
	    	 String s=getString(RGN_PARAM);
	    	 w.writeAttribute("value",s==null?"":s);

	    	 
	    	 w.writeEmptyElement("input");
	    	 w.writeAttribute("type", "submit");
	    	 w.writeAttribute("value","Submit");
	    	 w.writeEndElement();//div
	         w.writeEndElement();//form
	    	}
	    
	  
	    
	    private void doWork( XMLStreamWriter w,final GroupFile gf)
	    	throws XMLStreamException
	        {
	    	
	    	Position pos=parsePosition();
	        if(pos==null) return ;
	    	w.writeStartElement("div");

            w.writeStartElement("h2");
            w.writeCharacters(pos.chrom+":"+pos.pos);
            w.writeEndElement();
	    	
            w.writeStartElement("table");
            w.writeStartElement("thead");
            w.writeStartElement("tr");
    		for(String header:new String[]{"CHROM","POS","ID","REF","QUAL","Sample","Alleles",
    				"DP","GQ","File"})
	    		{
	    		w.writeStartElement("th");
	    		w.writeCharacters(header);
	    		w.writeEndElement();//td
	    		}

    		w.writeEndElement();//tr
    		w.writeEndElement();//thead
    		 w.writeStartElement("tbody");
    		 Set<String> samplesWithGenotypes=new HashSet<String>();
    		 Set<String> allSamples=new HashSet<String>();
        	for(VcfFile f:getVcfFiles(gf))
        		{
    			TabixReader tabixReader=null;
    			TabixReader.Iterator iter=null;
    			BlockCompressedInputStream bgzin=null;
    			VCFHeader header=null;
    		    AbstractVCFCodec vcfCodec = VCFUtils.createDefaultVCFCodec();
    		   LineIterator lineIterator=null;
    		   for(int i=0;i< 2;i++)
	    		   {
	        		try
	        			{
	        			if(i==0)
		        			{
			        		bgzin=new BlockCompressedInputStream(f.file);
	
		        			lineIterator=new LineIteratorImpl(new SynchronousLineReader(bgzin));
		        			header=(VCFHeader) vcfCodec.readActualHeader(lineIterator);
		        			allSamples.addAll(header.getGenotypeSamples());
		        			}
	        			else
	        				{
	        				tabixReader=new TabixReader(f.file.getPath());
	        				String line;
	        				
	        				int[] x = tabixReader.parseReg(pos.chrom+":"+pos.pos+"-"+(pos.pos));
	        				if(x!=null && x.length>2 && x[0]!=-1)
		        				{
		        				iter=tabixReader.query(x[0], x[1], x[2]);
		        				}
	        				else
	        					{
	        					}
	        				
	        				
	        				while(iter!=null && (line=iter.next())!=null)
	        					{
	        					VariantContext var=vcfCodec.decode(line);
	        					
	        					for(String sample:header.getSampleNamesInOrder())
	        						{
	        						final Genotype genotype=var.getGenotype(sample);
	        						if(genotype==null || !genotype.isCalled()) continue;
	        						if(!genotype.isAvailable()) continue;
	        						
	        						samplesWithGenotypes.add(sample);
	        						
	        						w.writeStartElement("tr");
	        						
	        						w.writeStartElement("td");
	        						w.writeCharacters(var.getContig());
	        						w.writeEndElement();

	        						
	        						w.writeStartElement("td");
	        						w.writeCharacters(String.valueOf(var.getStart()));
	        						w.writeEndElement();
	        						
	        						if(var.hasID())
	        							{
		        						w.writeStartElement("td");
		        						
		        							if( var.getID().matches("rs[0-9]+"))
		        								{
		        								w.writeStartElement("a");
		        								w.writeAttribute("href",
		        										"http://www.ncbi.nlm.nih.gov/snp/"+
		        										var.getID().substring(2)
		        										);
		        								w.writeCharacters(var.getID());
		        								w.writeEndElement();//a
		        								}
		        							else
		        								{
		        								w.writeCharacters(var.getID());
		        								}
		        						w.writeEndElement();//td
	        							}
	        						else
	        							{
	        							w.writeEmptyElement("td");
	        							}
	        						
	        						if(var.getReference()!=null)
		        						{
		        						w.writeStartElement("td");
		        						w.writeCharacters(var.getReference().getBaseString());
		        						w.writeEndElement();
		        						}
	        						else
	        							{
	        							w.writeEmptyElement("td");
	        							}
	        						
	        						if(var.hasLog10PError())
	        							{
		        						w.writeStartElement("td");
		        						w.writeCharacters(String.valueOf((int)var.getPhredScaledQual()));
		        						w.writeEndElement();
	        							}
	        						else
	        							{
	        							w.writeEmptyElement("td");
	        							}
	        						
	        						w.writeStartElement("td");
	        						w.writeCharacters(sample);
	        						w.writeEndElement();
	        						
	        						List<Allele> alleles=genotype.getAlleles();
	        						w.writeStartElement("td");
	        						
	        						w.writeStartElement("span");
	        						if(genotype.isHomRef())
	        							{
	        							w.writeAttribute("style", "color:green;");
	        							}
	        						else if(genotype.isHomVar())
	        							{
	        							w.writeAttribute("style", "color:red;");
	        							}
	        						else if(genotype.isHet())
	        							{
	        							w.writeAttribute("style", "color:blue;");
	        							}
	        						
	        						for(int j=0;j< alleles.size();++j)
	        							{
	        							if(j>0) w.writeCharacters(genotype.isPhased()?"|":"/");
	        							 w.writeCharacters(alleles.get(j).getBaseString());
	        							}
	        						w.writeEndElement();//span
	        						w.writeEndElement();
	        						
	        						
	        						if(genotype.hasDP())
	        							{
		        						w.writeStartElement("td");
		        						w.writeCharacters(String.valueOf(genotype.getDP()));
		        						w.writeEndElement();
	        							}
	        						else
	        							{
	        							w.writeEmptyElement("td");
	        							}
	        						
	        						if(genotype.hasGQ())
	        							{
		        						w.writeStartElement("td");
		        						w.writeCharacters(String.valueOf(genotype.getGQ()));
		        						w.writeEndElement();
	        							}
	        						else
	        							{
	        							w.writeEmptyElement("td");
	        							}
	        						w.writeStartElement("td");
	        						w.writeCharacters(f.file.getName());
	        						w.writeEndElement();

	        						
	        						w.writeEndElement();//tr
	        						w.flush();
	        						}
	        					}

	        				}
	        			}
	        		catch(Exception err)
	        			{
	        			w.writeComment("BOUM "+err);
	        			header=null;
	        			lastException=err;
	        			}
	        		finally
	        			{
	        			CloserUtil.close(lineIterator);
	        			CloserUtil.close(bgzin);
	        			CloserUtil.close(tabixReader);
	        			CloserUtil.close(iter);
	        			}
	        		if(i==0 && header==null) break;
	        		
	    		   }
    		   w.flush();
        		}
        	 w.writeEndElement();//tbody
        	 w.writeEndElement();//table
        	 
        	 allSamples.removeAll(samplesWithGenotypes);
        	 if(!allSamples.isEmpty())
        	 	{
                 w.writeStartElement("h3");
                 w.writeCharacters("Samples not found");
                 w.writeEndElement();
                 w.writeStartElement("ol");
                for(String sample:new TreeSet<String>(allSamples))
                	{
                    w.writeStartElement("li");
                    w.writeCharacters(sample);
                    w.writeEndElement();
                	}
                 w.writeEndElement();
        	 	}
        	 
	        writeHTMLException(w, this.lastException);

	    	w.writeEndElement();//div
	        }

	   
	    private void handleGroup(final GroupFile gf)
	        {
	    	setMimeHeaderPrinted(true);
	        System.out.print("Content-type: text/html;charset=utf-8\n");
	        System.out.println();
	        System.out.flush();
	       
	        Position rgn=parsePosition(getString(RGN_PARAM));
	        XMLStreamWriter w=null;
	    	try
	            {
	            XMLOutputFactory xof=XMLOutputFactory.newFactory();
	            w=xof.createXMLStreamWriter(System.out,"UTF-8");
	            w.writeStartElement("html");
	            w.writeStartElement("body");
	           
	            w.writeStartElement("h1");
	            w.writeCharacters(gf.getDesc());
	            w.writeEndElement();
	            printForm(w, gf);
	            
	            if(rgn!=null)
	            	{
	            	doWork(w,gf);
	            	}
	            else
	            	{
		            w.writeStartElement("h2");
		            w.writeCharacters("Files");
		            w.writeEndElement();
		            
		            w.writeStartElement("table");
		            
		            w.writeStartElement("tr");
            		
            		w.writeStartElement("th");
            		w.writeCharacters("File");
            		w.writeEndElement();//td

            		w.writeStartElement("td");
            		w.writeCharacters("Description");
            		w.writeEndElement();//td

            		w.writeEndElement();//tr
		            List<VcfFile> vcffiles=getVcfFiles(gf);
	            	for(VcfFile f:vcffiles)
	            		{
	            		w.writeStartElement("tr");
	            		
	            		w.writeStartElement("td");
	            		w.writeCharacters(f.file.getPath());
	            		w.writeEndElement();//td

	            		w.writeStartElement("td");
	            		w.writeCharacters(f.getDesc());
	            		w.writeEndElement();//td

	            		
	            		w.writeEndElement();//tr
	            		}
	            	 w.writeEndElement();//table
	            	}
	          
	           writeHTMLException(w, this.lastException);
	            writeHtmlFooter(w);
	            w.writeEndElement();//body
	            w.writeEndElement();//html
	            }
	        catch(Exception err)
	            {
	        	LOG.error(err);
	            }
	        finally
	            {
	            if(w!=null)
	                {
	                try {w.flush();} catch(XMLStreamException err){}
	                CloserUtil.close(w);
	                }
	            CloserUtil.close(System.out);
	            }
	    	
	    	

	        }
	   
	    @Override
	    protected void doCGI()
	        {
	        GroupFile gf=null;
	        String groupidstr=getString(GROUPID_PARAM);
	        if(groupidstr!=null)
	            {
	            try    {   
	                int groupid=Integer.parseInt(groupidstr);
	                List<GroupFile> gfs=getGroupFiles();
	                if(groupid>=0 && groupid<gfs.size())
	                    {
	                    gf=gfs.get(groupid);
	                    }
	                }
	            catch(Exception err)
	                {
	               
	                }
	            }
	        if(gf!=null)
	            {
	            handleGroup(gf);
	            return;
	            }
	        welcomePane();
	        }

	    /**
	     * @param args
	     */
	    public static void main(String[] args)
	        {
	        new VcfRegistryCGI().instanceMainWithExit(args);
	        }

}
