package com.github.lindenb.jvarkit.tools.ngsfiles;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.Date;
import java.util.HashSet;
import java.util.Set;

import javax.xml.stream.XMLEventFactory;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamResult;

import org.broadinstitute.variant.vcf.VCFHeader;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.illumina.FastQName;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.sleepycat.bind.tuple.StringBinding;
import com.sleepycat.je.Cursor;
import com.sleepycat.je.Database;
import com.sleepycat.je.DatabaseConfig;
import com.sleepycat.je.DatabaseEntry;
import com.sleepycat.je.Environment;
import com.sleepycat.je.EnvironmentConfig;
import com.sleepycat.je.LockMode;
import com.sleepycat.je.OperationStatus;
import com.sleepycat.je.Transaction;

public class NgsFilesScanner extends AbstractScanNgsFilesProgram
	{
	private static final String SUFFIXES[]=new String[]{".bam",".vcf",".vcf.gz"};
	private Transaction txn=null;
	private Database database=null;
	private FileFilter fileFilter=new FileFilter()
		{
		@Override
		public boolean accept(File f) {
			if(f==null) return false;
			if(!f.canRead()) return false;
			if(f.isDirectory()) return true;
			String filename=f.getName().toLowerCase();
			if(filename.contains("jeter")) return false;
			if(filename.contains("_delete_")) return false;
			for(String suff: SUFFIXES)
				if(filename.endsWith(suff)) return true;
			return false;
			}
		};
		
		private FileFilter fastqFilter=new FileFilter()
			{
			@Override
			public boolean accept(File f) {
				if(f==null) return false;
				if(!f.canRead()) return false;
				if(!f.isFile()) return true;
				String filename=f.getName().toLowerCase();
				if(filename.contains("jeter")) return false;
				if(filename.contains("_delete_")) return false;
				if(filename.endsWith(".fastq.gz")) return true;
				if(filename.endsWith(".fastq")) return true;
				if(filename.endsWith(".fq.gz")) return true;
				if(filename.endsWith(".fqz")) return true;
				return false;
				}
			};	
		
	static final String DATABASE_NAME="ngsfile.db";
	
    private NgsFilesScanner()
    	{
    	
    	}		
    @Override
    protected String getOnlineDocUrl() {
    	return "https://github.com/lindenb/jvarkit/wiki/NgsFilesScanner";
    	}
    
    @Override
    public String getProgramDescription() {
    	return "Build a persistent database of NGS file. Dump as XML. ";
    	}
    
    
    private void writeFile(XMLStreamWriter out,File f) throws XMLStreamException
    	{
    	out.writeAttribute("file", f.getAbsolutePath());
    	out.writeAttribute("filename", f.getName());
    	out.writeAttribute("size", String.valueOf(f.length()));
    	out.writeAttribute("modified", String.valueOf(new Date(f.lastModified())));
    	out.writeAttribute("timestamp", String.valueOf(f.lastModified()));
    	
    	}
    
    private void put(File f,String xml)
    	{
    	info("insert "+f);
		DatabaseEntry key=new DatabaseEntry();
		DatabaseEntry data=new DatabaseEntry();
		StringBinding.stringToEntry(f.getAbsolutePath(), key);
		StringBinding.stringToEntry(xml,data);
		this.database.put(this.txn, key, data);
    	}
    
    @Override
    protected void readFastq(File f) {
    	//do nothing
    	}
    
	@Override
    protected void readBam(File f)
    	{
    	if(!f.canRead()) return;
    	SAMFileReader r=null;
    	try {
    		StringWriter sw=new StringWriter();
    		XMLOutputFactory xof=XMLOutputFactory.newFactory();
    		XMLStreamWriter out=xof.createXMLStreamWriter(new StreamResult(sw));
    		
    		
    		out.writeStartElement("bam");
    		writeFile(out,f);
    		
			r=new SAMFileReader(f);
			r.setValidationStringency(ValidationStringency.LENIENT);
			SAMFileHeader h=r.getFileHeader();
			
			out.writeStartElement("samples");
			if(h!=null && h.getReadGroups()!=null)
				{
				Set<String> seen=new HashSet<String>();
				for(SAMReadGroupRecord rg: h.getReadGroups())
					{
					String sample=rg.getSample();
					if(sample==null || sample.isEmpty() || seen.contains(sample)) continue;
					
					seen.add(sample);
					out.writeStartElement("sample");
					out.writeCharacters(sample);
					out.writeEndElement();
					}
				}
			out.writeEndElement();
			
			out.writeEndElement();
			out.flush();
			out.close();
			sw.flush();
			
			put(f,sw.toString());
			} 
    	catch (Exception e)
    		{
    		warning(e, "Error in "+f);
			}
    	finally
    		{
    		if(r!=null) r.close();
    		}
    	}
   
   
    
   @Override
   protected void readVCF(File f)
		{
    	if(!f.canRead()) return;
    	debug("readVCF  "+f);
    	    	

    	VcfIterator r=null;
    	InputStream in=null;
    	try
    		{
    	
    		r=VCFUtils.createVcfIteratorFromFile(f);
        	VCFHeader header=r.getHeader();
        	
        	
    		StringWriter sw=new StringWriter();
    		XMLOutputFactory xof=XMLOutputFactory.newFactory();
    		XMLStreamWriter out=xof.createXMLStreamWriter(new StreamResult(sw));
    		out.writeStartElement("vcf");
    		writeFile(out,f);

    		out.writeStartElement("samples");
        	for(String sample:header.getSampleNamesInOrder())
	        	{
        		out.writeStartElement("sample");
				out.writeCharacters(sample);
				out.writeEndElement();
	    		}
        	out.writeEndElement();
        	
        	out.writeEndElement();
			out.flush();
			out.close();
			sw.flush();
			
			put(f,sw.toString());
    		}
    	catch(Exception err)
    		{
    		error(err,"Error in VCF "+f);
    		}
    	finally
    		{
    		CloserUtil.close(r);
    		CloserUtil.close(in);
    		}
    	
		}
   
   
   private void fastqDir(File dir,Counter<String> samples)
   	{
   	try {
   		StringWriter sw=new StringWriter();
   		XMLOutputFactory xof=XMLOutputFactory.newFactory();
   		XMLStreamWriter out=xof.createXMLStreamWriter(new StreamResult(sw));
   		
   		
   		out.writeStartElement("fastq-dir");
   		out.writeAttribute("directory", dir.getAbsolutePath());
		out.writeStartElement("samples");
		for(String sample:samples.keySet())
			{
			out.writeStartElement("sample");
			out.writeAttribute("size",String.valueOf(samples.count(sample)));
			out.writeCharacters(sample);
			out.writeEndElement();
			}
			
		out.writeEndElement();
		
		out.writeEndElement();
		out.flush();
		out.close();
		sw.flush();
		
		put(dir,sw.toString());
		} 
   	catch (Exception e)
   		{
   		warning(e, "Error in "+dir);
		}
   	finally
   		{
   		}
   	
   }
   
   private static final  Counter<String> EMPTY_COUNTER=new Counter<String>();
   
   private Counter<String> getFastqSampleInDirectory(File f)
	   	{
	   if(f==null || !f.isDirectory() || !f.exists() || !f.canRead()) return EMPTY_COUNTER;
	    File array[]=f.listFiles(this.fastqFilter);
	   	if(array==null) return EMPTY_COUNTER;
	   	Counter<String> fastqSamples=new Counter<String>();
	   	for(File f2:array)
	   		{
	   		FastQName fqName=FastQName.parse(f2);
			if(fqName.isValid() && fqName.getSample()!=null)
				{
				fastqSamples.incr(fqName.getSample(),f2.length());
				}
	   		}
	   	return fastqSamples;
	   	}
	
    private long lastDirTimeMillis=System.currentTimeMillis();
    private void recursive(File f) throws IOException
    	{
    	if(f==null || !f.exists() || !f.canRead()) return;
    	
    	if(f.getName().startsWith(".")) return;
    	
    	if(f.isDirectory() && f.canRead())
    		{
    		if(f.getName().toLowerCase().equals("tmp")) return;
    		if(f.getName().toLowerCase().equals("jeter")) return;
    		
    		
    		if(f.getName().equals("Intensities") && f.getParentFile()!=null &&
    				f.getParentFile().getName().equals("Data"))
    			{
    			info("Skipping "+f);
    			return;
    			}
    		if(f.getName().startsWith("L") && f.getParentFile()!=null &&
    				(f.getParentFile().getName().equals("Thumbnail_Images") || 
    				 f.getParentFile().getName().equals("Processed")
    				))
    			{
    			info("Skipping "+f);
    			return;
    			}
    		
    		long now=System.currentTimeMillis();
    		if(now-lastDirTimeMillis > 30*1000)
	    		{
    			info("In "+f);
	    		this.lastDirTimeMillis=now;
	    		}
    		
    		
    		File array[]=f.listFiles(this.fileFilter);
    		if(array!=null)
    			{
    			for(File f2:array) recursive(f2);    		
    			}
    		
    		
    		Counter<String> fastqSamples=getFastqSampleInDirectory(f);
    		if(!fastqSamples.isEmpty())
    			{
    			fastqDir(f, fastqSamples);
    			}
    		}
    	else if(f.isFile() && this.fileFilter.accept(f))
    		{
    		analyze(f);
    		}
    	}
    
    private boolean isEntryShouldBeDeleted(File f)
    	{
    	if(!f.exists()) return true;
    	if(f.isDirectory() && getFastqSampleInDirectory(f).isEmpty()) return true;
    	return false;
    	}
    
    @Override
    public void printOptions(PrintStream out) {
    	out.println("-B (dir) "+getMessageBundle("berkeley.db.home"));
    	out.println("-D dump as XML to stdout and exit");
    	super.printOptions(out);
    	}

	@Override
	public int doWork(String[] args)
		{
		boolean dump=false;
		EnvironmentConfig envCfg=new EnvironmentConfig();
		Environment env=null;
		File bdbHome=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"B:D"))!=-1)
			{
			switch(c)
				{
				case 'D': dump=true;break;
				case 'B': bdbHome=new File(opt.getOptArg());break;
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
		
		
		if(bdbHome==null)
			{
			error("BDB home undefined");
			return -1;
			}
		if(!bdbHome.exists() || !bdbHome.isDirectory())
			{
			error("BDB doesn't exist or is not a directory");
			return -1;
			}
		
		File root=new File("/");
		Cursor cursor=null;
		try
			{
			
			if(!dump && opt.getOptInd()+1==args.length)
				{
				root=new File(args[opt.getOptInd()]);
				}
			else if(!(dump && opt.getOptInd()==args.length))
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			
			envCfg.setAllowCreate(!dump);
			envCfg.setReadOnly(dump);
			envCfg.setTransactional(false);
			
			
			info("Opening env "+bdbHome);
			env=new Environment(bdbHome, envCfg);
			
			
			//TransactionConfig txnCfg=new TransactionConfig();
			//this.txn=env.beginTransaction(null, txnCfg);

			
			info("Opening database "+DATABASE_NAME);
			DatabaseConfig cfg=new DatabaseConfig();
			cfg.setAllowCreate(!dump);
			cfg.setReadOnly(dump);
			cfg.setTransactional(false);
			this.database=env.openDatabase(this.txn,DATABASE_NAME, cfg);
			
			DatabaseEntry key=new DatabaseEntry();
			DatabaseEntry data=new DatabaseEntry();

			if(dump)
				{
				XMLEventFactory xef=XMLEventFactory.newFactory();
	    		XMLOutputFactory xof=XMLOutputFactory.newFactory();
	    		XMLInputFactory xif=XMLInputFactory.newFactory();
	    		
	    		XMLEventWriter out=xof.createXMLEventWriter(System.out,"UTF-8");
	    		out.add(xef.createStartDocument("UTF-8", "1.0"));
	    		out.add(xef.createStartElement("", "","ngs-files"));
	    		cursor=this.database.openCursor(this.txn, null);
	    		while(cursor.getNext(key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS)
	    			{
	    			File file=new File(StringBinding.entryToString(key));
	    			if(isEntryShouldBeDeleted(file))
	    				{
	    				info("deleting entry for "+file);
	    				//cursor.delete();//no env is read only
	    				continue;
	    				}
	    			StringReader sr=new StringReader(StringBinding.entryToString(data));
	    			XMLEventReader xr=xif.createXMLEventReader(sr);
	    			while(xr.hasNext())
	    				{
	    				XMLEvent evt=xr.nextEvent();
	    				if(evt.isStartDocument() ) continue;
	    				if(evt.isEndDocument()) break;
	    				out.add(evt);
	    				}
	    			xr.close();
	    			sr.close();
	    			}
	    		cursor.close();
	    		
	    		out.add(xef.createEndElement("","","ngs-files"));
	    		out.add(xef.createEndDocument());
	    		out.flush();
	    		out.close();
	    		System.out.flush();
	    		System.out.close();
				}
			else
				{
				recursive(root);
				
				//final cleanup
				cursor=this.database.openCursor(this.txn, null);
				while(cursor.getNext(key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS)
	    			{
	    			File file=new File(StringBinding.entryToString(key));
	    			
	    			if(isEntryShouldBeDeleted(file))
	    				{
	    				info("deleting entry for "+file);
	    				cursor.delete();
	    				}
	    			}
				cursor.close();
				}
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			
			if(cursor!=null)  try { cursor.close();} catch(Exception err){}
			if(this.txn!=null)try { this.txn.commit();} catch(Exception err){}
			if(this.database!=null) try { this.database.close();} catch(Exception err){}
			if(env!=null) try { env.close();} catch(Exception err){}
			}
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		/* https://bugs.openjdk.java.net/browse/JDK-8028111 */
		System.setProperty("jdk.xml.entityExpansionLimit","0");
		new NgsFilesScanner().instanceMainWithExit(args);

	}

}
