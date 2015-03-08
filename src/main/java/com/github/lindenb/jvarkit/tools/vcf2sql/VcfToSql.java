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
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;

public class VcfToSql extends AbstractCommandLineProgram
	{
    
    private static long ID_GENERATOR=0L;
    private static final String NS="http://github.com/lindenb/jvarkit/";
    
    
    
    private abstract class AbstractComponent
    	{
    	String name;
    	String rdfsLabel=null;
    	String rdfsComment=null;
    	String d2rquriPattern=null;
    	
    	AbstractComponent(String name)
			{
			this.name=name;
			}
    	
    	AbstractComponent label(String s) {this.rdfsLabel=s;return this;}
    	AbstractComponent comment(String s) {this.rdfsComment=s;return this;}
    	AbstractComponent uriPattern(String s) {this.d2rquriPattern=s;return this;}
    	
    	public String getLabel() { return this.rdfsLabel==null?getName():this.rdfsLabel;}
    	public String getComment() { return this.rdfsComment==null?getLabel():this.rdfsComment;}
    	public String getUriPattern() { return d2rquriPattern;}
    	public String getD2RQName()
    		{
    		return this.getName();
    		}

    	public String getName()
    		{
			return name;
			}
    	public String getAntiquote()
    		{
			return "`"+getName()+"`";
			}
    	
    	
   	
    	public abstract void createD3RQMapping(PrintWriter pw);
    	}
    
    private class ColumnBuilder
    	{
    	private Table _foreignTable=null;
    	private String _comment=null;
    	private String _label=null;
    	private boolean _pkey=false;
    	private String _name="";
    	private Class _class=String.class;
    	ColumnBuilder foreignKey(Table v) { this._foreignTable=v; return this;}
    	ColumnBuilder primaryKey() { this._pkey=true; return this;}
    	ColumnBuilder name(String v) { this._name=v; return this;}
    	ColumnBuilder comment(String v) { this._comment=v; return this;}
    	ColumnBuilder label(String v) { this._label=v; return this;}
    	ColumnBuilder type(Class v) { this._class=v; return this;}
    	
    	Column make()
    			{
    			Column c=null;
    			if(_foreignTable!=null)
    				{
    				c = new ForeignKey(this._foreignTable);
    				}
    			else if(_pkey)
    				{
    				c = new PrimaryKey();
    				}
    			else if(_class==Double.class)
    				{
    				c = new DoubleColumn(this._name);
    				}
    			else if(_class==Integer.class)
    				{
    				c = new IntegerColumn(this._name);
    				}
    			else  if(_class==String.class)
    				{
    				c = new StringColumn(this._name);
    				}
    			else
    				{
    				throw new IllegalArgumentException();
    				}
    			
    			c.rdfsLabel = (this._label==null?c.getName():this._label);
    			c.rdfsComment = (this._comment==null?c.rdfsLabel:this._label);
    			
    			return c;
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
    			table.out.print("\\N");
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
    	
    	public abstract String mysqlCreateStatement();
    	
    	public String getD2RQName()
    		{
    		return table.getD2RQName()+"_"+this.getName();
    		}

    	
    	@Override
    	public void createD3RQMapping(PrintWriter pw)
    		{
    		pw.println("map:"+getD2RQName()+" a d2rq:PropertyBridge;");
    		pw.println("	d2rq:belongsToClassMap map:"+table.getName()+";");
    		pw.println("	d2rq:column \""+table.getName()+"."+this.getName()+"\";");
    		pw.println("	d2rq:property	"+this.getName()+";");
    		pw.println("	d2rq:propertyDefinitionLabel \""+getLabel()+"\"@en ;");
    		pw.println("	d2rq:propertyDefinitionComment \""+getComment()+"\"@en ;");
    		pw.println("	.");
    		}

    	
    	}
    
    private class LongColumn
	extends Column
		{
		LongColumn(String name)
			{
			super(name);
			}
		@Override
		public String mysqlCreateStatement()
			{
			return getAntiquote()+" INT "+(nilleable?"":" NOT ")+"NULL";
			}
		}

    
    private class ForeignKey
	extends LongColumn
		{
		Table references;
    	public ForeignKey(Table t)
    		{
			super(t.getName()+"_id");
			this.references=t;
    		}
    	@Override
    	public String mysqlCreateStatement()
    		{
    		return getAntiquote()+" INT NOT NULL, "+
    			"FOREIGN KEY ("+getAntiquote()+")  " +
    			"REFERENCES "+references.getAntiquote()+"(id)"
    			; 
    		}

    	@Override
       	public void createD3RQMapping(PrintWriter pw)
       		{
    		pw.println("map:"+getD2RQName()+"_ref a d2rq:PropertyBridge;");
    		pw.println("	d2rq:belongsToClassMap map:"+this.table.getD2RQName()+";");
    		pw.println("	d2rq:property vocab:"+this.references.getD2RQName()+"_id;");
    		pw.println("	d2rq:refersToClassMap map:"+this.references.getD2RQName()+";");
    		pw.println("	d2rq:join \""+this.table.getName()+"."+this.getName()+" => "+ references.getName() +".id\";");
    		pw.println("	.");
       		}

		}
    
    private class PrimaryKey
   	extends LongColumn
   		{
   		String d2rq_pattern=null;
   		
       	public PrimaryKey()
       		{
   			super("id");
       		}
       	
       	boolean isPrimaryKey()
       		{
       		return true;
       		}
       	@Override
    	public String mysqlCreateStatement()
    		{
    		return getAntiquote()+" INT NOT NULL, "+
    			"PRIMARY KEY ("+getAntiquote()+")" 
    			; 
    		}
       	@Override
       	public void createD3RQMapping(PrintWriter pw)
       		{
       		pw.print("map:"+getD2RQName()+" a d2rq:PropertyBridge;");
       		pw.print("	d2rq:belongsToClassMap map:"+ table.getD2RQName() +";");
       		pw.print("	d2rq:property rdfs:label;");
       		if(d2rq_pattern==null)
       			{
       			pw.print("	d2rq:pattern \""+ table.getName()+"#@@"+ table.getName()+".id@@\";");
       			}
       		else
       			{
       			pw.print("	d2rq:pattern \""+ this.d2rq_pattern +"\";");
       			}
       		pw.print(".");

       		}
   		}
    
    private class IntegerColumn
	extends Column
		{
    	IntegerColumn(String name)
			{
			super(name);
			}
		@Override
		public String mysqlCreateStatement()
			{
			return getAntiquote()+" INT "+(nilleable?"":" NOT ")+"NULL";
			}

		}
    
    private class DoubleColumn
	extends Column
		{
    	DoubleColumn(String name)
			{
			super(name);
			}
		@Override
		public String mysqlCreateStatement()
			{
			return getAntiquote()+" DOUBLE "+(nilleable?"":" NOT ")+"NULL";
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
				table.out.print("\\N");
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
		@Override
		public String mysqlCreateStatement()
			{
			return getAntiquote()+" VARCHAR("+(maxLength+1)+") "+
					(nilleable?"":" NOT ")+"NULL";
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
    	String owlClass="rdfs:Class";
    	
    	Table(String name,Column...columns)
			{
			super(name);
			this.columns.addAll(Arrays.asList(columns));
			for(Column c:columns)
				{
				c.table=this;
				this.name2column.put(c.getName(), c);
				}
			super.d2rquriPattern=NS+name+"/@@"+name+".id@@";
			}
    	
    	Table owlClass(String owlClass) { this.owlClass=owlClass; return this; }
    	
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
    		this.last_inserted_id=Number.class.cast(row[0]).longValue();
    		
    		for(int i=0;i < this.columns.size();++i)
    			{
    			if(i>0 ) this.out.print(this.fieldDelimiter);
    			this.columns.get(i).print(row[i]);
    			}
    		this.out.println();
    		this.nRows++;
    		return this.last_inserted_id;
    		}
    	
    	public String mysqlCreateStatement()
    		{
    		StringBuilder sw=new StringBuilder();
    		sw.append("CREATE TABLE IF NOT EXISTS `"+getName()+"` (" );
    		if(!columns.get(0).isPrimaryKey())
    			{
    			sw.append("`id` INT NOT NULL AUTO_INCREMENT");
    			}
    		
    		for(int i=0;i< columns.size();++i)
    			{
    			if(i>0 || !columns.get(0).isPrimaryKey())
    				{
    				sw.append(",");
    				}
    			sw.append(columns.get(i).mysqlCreateStatement());
    			}
    		
    		if(!columns.get(0).isPrimaryKey())
    			{
    			sw.append(", PRIMARY KEY (`id`)");
    			}
    		
    		
    		sw.append(") DEFAULT CHARSET=utf8 ;");
    		return sw.toString();
    		}
    	/* http://plindenbaum.blogspot.fr/2010/02/mapping-rdbms-to-rdf-with-d2rq-yet.html */
    	@Override
    	public void createD3RQMapping(PrintWriter pw)
    		{
    		pw.println("\n# Table "+getName());
    		pw.println("map:"+getD2RQName()+" a d2rq:ClassMap;");
    		pw.println("	d2rq:uriPattern	\""+this.d2rquriPattern +"\" ;");
    		pw.println("	d2rq:class	"+this.owlClass+" ;");
    		pw.println("	d2rq:classDefinitionLabel \""+getLabel()+"\"@en ;");
    		pw.println("	d2rq:classDefinitionComment \""+getComment()+"\"@en ;");
    		pw.println("	d2rq:dataStorage map:Database1;");    		
    		pw.println("	.");
    		
       		pw.print("map:"+this.getD2RQName()+"__label a d2rq:PropertyBridge;");
       		pw.print("	d2rq:belongsToClassMap map:"+ this.getD2RQName() +";");
       		pw.print("	d2rq:property rdfs:label;");
       		pw.print("	d2rq:pattern \""+ this.getName()+"#@@"+ this.getName()+".id@@\";");
       		pw.print(".");
       		
       
    		for(Column c:this.columns)
    			{
    			c.createD3RQMapping(pw);
    			}
    		}
    	}
    
    
    private Table vcfFileTable = new Table("vcffile",
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().name("file").make()
    		);
    private Table sampleTable = new Table("sample",
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(vcfFileTable).make(),
    		new ColumnBuilder().name("name").make()
    		);
    
    private Table filterTable = new Table("filter",
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(vcfFileTable).make(),
    		new ColumnBuilder().name("name").make(),
    		new ColumnBuilder().name("description").make()
    		);

    private Table chromosomeTable = new Table("chromosome",
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(vcfFileTable).make(),
    		new ColumnBuilder().name("name").make(),
    		new ColumnBuilder().name("chromLength").type(Integer.class).make()
    		);
    
    private Table variantTable = new Table("variant",
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(vcfFileTable).make(),
    		new ColumnBuilder().name("index_in_file").type(Integer.class).make(),
    		new ColumnBuilder().foreignKey(chromosomeTable).make(),
    		new ColumnBuilder().name("pos").type(Integer.class).make(),
    		new ColumnBuilder().name("rsid").make(),
    		new ColumnBuilder().name("ref").make(),
    		new ColumnBuilder().name("pos").type(Double.class).make()
    		); 
    
    private Table variant2altTable = new Table("variant2alt",
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(variantTable).make(),
    		new ColumnBuilder().name("alt").make()
    		); 
    
    private Table variant2filters = new Table("variant2filter",
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(variantTable).make(),
    		new ColumnBuilder().foreignKey(filterTable).make()
    		); 

    private Table vepPrediction = new Table("vepPrediction",
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(variantTable).make(),
    		new ColumnBuilder().name("ensGene").make(),
    		new ColumnBuilder().name("ensTranscript").make(),
    		new ColumnBuilder().name("ensProtein").make(),
    		new ColumnBuilder().name("geneSymbol").make()
    		); 
    private Table vepPrediction2so = new Table("vepPrediction2so",
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(vepPrediction).make(),
    		new ColumnBuilder().name("acn").make()
    		); 

    
    private Table genotypeTable = new Table("genotype",
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(variantTable).make(),
    		new ColumnBuilder().foreignKey(sampleTable).make(),
    		new ColumnBuilder().name("a1").make(),
    		new ColumnBuilder().name("a2").make(),
    		new ColumnBuilder().name("dp").type(Integer.class).make(),
    		new ColumnBuilder().name("gq").type(Double.class).make()
    		); 

    private Table all_tables[]=new Table[]{
    	   vcfFileTable,sampleTable,filterTable,chromosomeTable,
    	   variantTable,variant2altTable,variant2filters,
    	   vepPrediction,
    	   vepPrediction2so,
    	   genotypeTable
    	};
    
    
    public VcfToSql()
    	{
    	
    	}
    
    private static String emitN3(String subject,String...nodes)
		{
		StringBuilder sb=new StringBuilder(subject);
		for(int i=0;i< nodes.length;i+=2)
			{
			sb.append(nodes[i]).append(" ").append(nodes[i+1]);
			sb.append("\n\t");
			sb.append(i+2==nodes.length?".":";");
			}
		
		return sb.toString();
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -o (out.zip) filename out.");
		super.printOptions(out);
		}

	private long nextKey()
		{
		return ++ID_GENERATOR;
		}
	
	private void read(File filename)
		throws IOException
		{

		/* insert this sample */
		final long vcffile_id = this.vcfFileTable.insert(nextKey(),filename);
		Map<String,Long> sample2sampleid = new HashMap<String,Long>();
		Map<String,Long> filter2filterid = new HashMap<String,Long>();
		
		VcfIterator r=VCFUtils.createVcfIteratorFromFile(filename);
		VCFHeader header=r.getHeader();
		
		/* parse samples */
		for(String sampleName:header.getSampleNamesInOrder())
			{
			long sampleid = this.sampleTable.insert(nextKey(),vcffile_id,sampleName);
			sample2sampleid.put(sampleName, sampleid);
			}
		
		/* parse filters */
		for(VCFFilterHeaderLine filter:header.getFilterLines())
			{
			long filterid = this.filterTable.insert(
					nextKey(),
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
					(long)ssr.getSequenceIndex(),
					vcffile_id,
					ssr.getSequenceName(),
					ssr.getSequenceLength()
					);
			}
		
		VepPredictionParser vepPredictionParser=new VepPredictionParser(header);
		
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
		int nVariants=0;
		while(r.hasNext())
			{
			VariantContext var= progress.watch(r.next());
			++nVariants;
			
			/* insert variant */
			final long variant_id = this.variantTable.insert(
				nextKey(),
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
				this.variant2altTable.insert(
					nextKey(),
					variant_id,
					alt.getBaseString()
					);
				}
			
			/* insert filters */
			for(String filter:var.getFilters())
				{
				this.variant2filters.insert(
					nextKey(),
					variant_id,
					filter2filterid.get(filter)
					);
				}
			
			for(VepPrediction pred: vepPredictionParser.getPredictions(var))
				{
				long pred_id= vepPrediction.insert(
						nextKey(),
						variant_id,
						pred.getEnsemblGene(),
						pred.getEnsemblTranscript(),
						pred.getEnsemblProtein(),
						pred.getGene()
						);
						
				for(SequenceOntologyTree.Term t: pred.getSOTerms())
					{
					vepPrediction2so.insert(
						nextKey(),
						pred_id,
						t.getAcn()
						);
					}
				}
			
			/* insert genotypes */
			for(String sampleName: sample2sampleid.keySet())
				{
				Genotype g= var.getGenotype(sampleName);
				
				genotypeTable.insert(
						nextKey(),
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
	

	@Override
	public int doWork(String[] args)
		{
		final String zipDir="vcf2sql.output/";
		File zipFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:"))!=-1)
			{
			switch(c)
				{
				case 'o':
					String filename=opt.getOptArg();
					if(!filename.endsWith(".zip"))
						{
						error(filename +" must end with '.zip'");
						return -1 ;
						}
					zipFile=new File(filename);
					if(zipFile.getParentFile()!=null)
							this.addTmpDirectory(zipFile.getParentFile());
					break;
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
			if(opt.getOptInd()+1!=args.length)
				{
				info("Illegal number of arguments");
				return -1;
				}
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
				ZipEntry zipEntry=new ZipEntry(zipDir+t.getName()+".tsv");
				zout.putNextEntry(zipEntry);
				InputStream tmpIn = new FileInputStream(t.tmpFile);
				IOUtils.copyTo(tmpIn, zout);
				tmpIn.close();
				zout.closeEntry();
				t.tmpFile.delete();
				}
			
			/** schema */
			zout.putNextEntry(new ZipEntry(zipDir+"truncate.mysql"));
			PrintWriter pw=new PrintWriter(zout);
			for(Table t:this.all_tables)
				{
				pw.println("truncate TABLE "+t.getAntiquote()+";");
				}
			pw.flush();
			zout.closeEntry();
			
			zout.putNextEntry(new ZipEntry(zipDir+"drop.mysql"));
			pw=new PrintWriter(zout);
			for(int i=this.all_tables.length-1;i>=0;--i)
				{
				pw.println("DROP TABLE IF EXISTS "+all_tables[i].getAntiquote()+";");
				}
			pw.flush();
			zout.closeEntry();
			
			zout.putNextEntry(new ZipEntry(zipDir+"load.mysql"));
			pw=new PrintWriter(zout);
			for(Table t:this.all_tables)
				{
				pw.println("LOAD DATA LOCAL INFILE '"+t.getName()+".tsv' " +
						" INTO TABLE "+t.getAntiquote()+
						" FIELDS TERMINATED BY '\\t' ;");
				}
			pw.flush();
			zout.closeEntry();
			
			zout.putNextEntry(new ZipEntry(zipDir+"create.mysql"));
			pw=new PrintWriter(zout);
			for(Table t:this.all_tables)
				{
				pw.println(t.mysqlCreateStatement());
				}
			pw.flush();
			zout.closeEntry();

			zout.putNextEntry(new ZipEntry(zipDir+"mapping.n3"));
			pw=new PrintWriter(zout);
			pw.println("# D2RQ Namespace");  
			pw.println("@prefix d2rq:        <http://www.wiwiss.fu-berlin.de/suhl/bizer/D2RQ/0.1#> .");  
			pw.println("# Namespace of the ontology");  
			pw.println("@prefix : <http://annotation.semanticweb.org/iswc/iswc.daml#> .");  
			pw.println("# Namespace of the mapping file; does not appear in mapped data");  
			pw.println("@prefix map: <file:///Users/d2r/example.ttl#> .");  
			pw.println("# Other namespaces"); 
			pw.println("@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> ."); 
			pw.println("@prefix xsd: <http://www.w3.org/2001/XMLSchema#> . "); 
			pw.println("map:Database1 a d2rq:Database;"); 
			pw.println("	d2rq:jdbcDSN \"JDBC_URI\";"); 
			pw.println("	d2rq:jdbcDriver \"com.mysql.jdbc.Driver\";"); 
			pw.println("	d2rq:username \"JDBC_USER\";"); 
			pw.println("	d2rq:password \"JDBC_PASSWORD\"."); 
			
			for(Table t:this.all_tables)
				{
				t.createD3RQMapping(pw);
				}
			pw.flush();
			zout.closeEntry();

			
			zout.finish();
			zout.flush();
			zout.close();
			
			
			
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
    


	public static void main(String[] args)
		{
		new VcfToSql().instanceMainWithExit(args);
		}
	}
