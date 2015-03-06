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
package com.github.lindenb.jvarkit.tools.vcf2sql;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfToSql extends AbstractCommandLineProgram
	{
    
    private static long ID_GENERATOR=0L;
    
    
    private abstract class AbstractComponent
    	{
    	String name;
    	
    	AbstractComponent(String name)
			{
			this.name=name;
			}
    	
    	public String getName()
    		{
			return name;
			}
    	
    	}
    
    private abstract class Column
    	extends AbstractComponent
    	{
    	boolean nilleable=false;
    	Table table;
    	Column(String name)
			{
			super(name);
			}
    	void print(Object o)
    		{
    		if(o==null)
    			{
    			table.out.print("NULL");
    			nilleable=true;
    			}
    		else
    			{
    			table.out.print(String.valueOf(o));
    			}
    		}
    	boolean isPrimaryKey()
    		{
    		return false;
    		}
    	}
    
    private class ForeignKey
	extends Column
		{
    	public ForeignKey(Table t)
    		{
			super(t.getName()+"_id");
    		}
		}
    
    private class PrimaryKey
   	extends Column
   		{
       	public PrimaryKey()
       		{
   			super("id");
       		}
       	boolean isPrimaryKey()
       		{
       		return true;
       		}
   		}
    
    private class IntegerColumn
	extends Column
		{
    	IntegerColumn(String name)
			{
			super(name);
			}
		}
    
    private class DoubleColumn
	extends Column
		{
    	DoubleColumn(String name)
			{
			super(name);
			}
		}

    
    private class StringColumn
    	extends Column
    	{
    	int maxLength=0;
    	StringColumn(String name)
			{
			super(name);
			}
    	void print(Object o)
			{
    		if(o==null)
				{
				table.out.print("NULL");
				nilleable=true;
				}
			else
				{
				String s=String.valueOf(o);
				maxLength=Math.max(maxLength, s.length());
				for(int i=0;i< s.length();++i)
					{
					if(s.charAt(i)=='\'') table.out.print("\'");
					table.out.print(s.charAt(i));
					}
				}
			}
    	}
    
    private class Table
    	extends AbstractComponent
    	{
    	char fieldDelimiter='\t';
    	long nRows=0L;
    	File tmpFile=null;
    	long last_inserted_id=0L;
    	PrintWriter out=null;
    	List<Column> columns=new ArrayList<>();
    	Map<String,Column> name2column=new HashMap<String,Column>();
    	
    	Table(String name,Column...columns)
			{
			super(name);
			this.columns.addAll(Arrays.asList(columns));
			for(Column c:columns)
				{
				c.table=this;
				this.name2column.put(c.getName(), c);
				}
			}
    	
    	void open() throws IOException
    		{
    		if(out!=null) throw new IllegalStateException(); 
    		this.tmpFile = File.createTempFile(
    				"_vcf2sql."+this.getName(),".txt",
    				VcfToSql.this.getTmpDirectories().get(0)
    				); 
    		this.out = new PrintWriter(this.tmpFile);
    		}
    	
    	long insert(Object...row)
    		{
    		this.last_inserted_id=++ID_GENERATOR;
    		if(row.length!=this.columns.size())
    			{
    			throw new IllegalStateException();
    			}
    		this.out.print(last_inserted_id);
    		for(int i=0;i<this.columns.size();++i)
    			{
    			this.out.print(this.fieldDelimiter);
    			this.columns.get(i).print(row[i]);
    			}
    		this.out.println();
    		
    		++nRows;
    		return last_inserted_id;
    		}
    	}
    
    
    private Table vcfFileTable = new Table("vcffile",
    		new StringColumn("file")
    		);
    private Table sampleTable = new Table("sample",
    		new ForeignKey(vcfFileTable),
    		new StringColumn("name")
    		);
    
    private Table filterTable = new Table("filter",
    		new ForeignKey(vcfFileTable),
    		new StringColumn("name"),
    		new StringColumn("description")
    		);

    private Table chromosomeTable = new Table("chromosome",
    		new PrimaryKey(),
    		new ForeignKey(vcfFileTable),
    		new StringColumn("name"),
    		new IntegerColumn("chromLength")
    		);
    
    private Table variantTable = new Table("variant",
    		new ForeignKey(vcfFileTable),
    		new IntegerColumn("index_in_file"),
    		new ForeignKey(chromosomeTable),
    		new IntegerColumn("pos"),
    		new StringColumn("rsid"),
    		new StringColumn("ref"),
    		new DoubleColumn("qual")
    		); 
    
    private Table variant2altTable = new Table("variant2alt",
    		new ForeignKey(variantTable),
    		new StringColumn("alt")
    		); 
    
    private Table variant2filters = new Table("variant2filter",
    		new ForeignKey(variantTable),
    		new ForeignKey(filterTable)
    		); 

    private Table genotypeTable = new Table("genotype",
    		new ForeignKey(variantTable),
    		new ForeignKey(sampleTable),
    		new StringColumn("a1"),
    		new StringColumn("a2"),
    		new IntegerColumn("dp"),
    		new DoubleColumn("gq")
    		); 

    private Table all_tables[]=new Table[]{
    	   vcfFileTable,sampleTable,filterTable,chromosomeTable,
    	   variantTable,variant2altTable,variant2filters,genotypeTable
    	};
    
    private PrintWriter out=new PrintWriter(System.out);
    
    public VcfToSql()
    	{
    	
    	}
    
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -o (out.zip) filename out.");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		File zipFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:"))!=-1)
			{
			switch(c)
				{
				case 'o': zipFile=new File(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		
		try
			{
			if(opt.getOptInd()+1==args.length)
				{
				File filename=new File(args[opt.getOptInd()]);
				for(Table t:this.all_tables)
					{
					t.open();
					}
				read(filename);
				FileOutputStream fout=new FileOutputStream(zipFile);
				ZipOutputStream zout = new ZipOutputStream(fout);
				for(Table t:this.all_tables)
					{
					t.out.flush();
					t.out.close();
					ZipEntry zipEntry=new ZipEntry(t.getName()+".tsv");
					zout.putNextEntry(zipEntry);
					InputStream tmpIn = new FileInputStream(t.tmpFile);
					IOUtils.copyTo(tmpIn, zout);
					tmpIn.close();
					zout.closeEntry();
					t.tmpFile.delete();
					}
				zout.finish();
				zout.flush();
				zout.close();
				}
			else
				{
				info("Illegal number of arguments");
				return -1;
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
			
			}
		}
    
	
	
	private void read(File filename)
		throws IOException
		{

		/* insert this sample */
		final long vcffile_id = this.vcfFileTable.insert(filename);
		Map<String,Long> sample2sampleid = new HashMap<String,Long>();
		Map<String,Long> filter2filterid = new HashMap<String,Long>();
		
		VcfIterator r=VCFUtils.createVcfIteratorFromFile(filename);
		VCFHeader header=r.getHeader();
		
		/* parse samples */
		for(String sampleName:header.getSampleNamesInOrder())
			{
			long sampleid = this.sampleTable.insert(vcffile_id,sampleName);
			sample2sampleid.put(sampleName, sampleid);
			}
		
		/* parse filters */
		for(VCFFilterHeaderLine filter:header.getFilterLines())
			{
			long filterid = this.filterTable.insert(
					vcffile_id,
					filter.getID(),
					filter.getValue()
					);
			filter2filterid.put(filter.getID(), filterid);
			}

		
		final SAMSequenceDictionary dict= header.getSequenceDictionary();
		/* parse sequence dict */
		for(SAMSequenceRecord ssr: dict.getSequences())
			{
			this.chromosomeTable.insert(
					ssr.getSequenceIndex(),
					vcffile_id,
					ssr.getSequenceName(),
					ssr.getSequenceLength()
					);
			}
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
		int nVariants=0;
		while(r.hasNext())
			{
			VariantContext var= progress.watch(r.next());
			++nVariants;
			
			/* insert variant */
			final long variant_id = this.variantTable.insert(
				vcffile_id,
				nVariants,
				dict.getSequence(var.getChr()).getSequenceIndex(),
				var.getStart(),
				(var.hasID()?var.getID():null),
				var.getReference().getBaseString(),
				(var.hasLog10PError()?var.getPhredScaledQual():null)
				);
			
			/* insert alternate alleles */
			for(Allele alt: var.getAlternateAlleles())
				{
				this.variant2altTable.insert(variant_id,
					alt.getBaseString()
					);
				}
			
			/* insert filters */
			for(String filter:var.getFilters())
				{
				this.variant2filters.insert(
					variant_id,
					filter2filterid.get(filter)
					);
				}
			
			/* insert genotypes */
			for(String sampleName: sample2sampleid.keySet())
				{
				Genotype g= var.getGenotype(sampleName);
				
				genotypeTable.insert(
						variant_id,
						sample2sampleid.get(sampleName),
						g.isCalled()?g.getAllele(0).getBaseString():null,
						g.isCalled()?g.getAllele(1).getBaseString():null,
						g.hasDP()?g.getDP():null,
						g.hasGQ()?g.getGQ():null	
						);
				}
			
			}
		r.close();
		}
	


	public static void main(String[] args)
		{
		new VcfToSql().instanceMainWithExit(args);
		}
	}
