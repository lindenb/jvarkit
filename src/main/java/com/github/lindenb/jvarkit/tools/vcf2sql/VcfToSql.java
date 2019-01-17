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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.vcf2sql;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

/**
BEGIN_DOC

## Examples

```bash
java -jar dist/vcf2sql.jar  file.vcf | mysql -u user -p -D vcf_db 
```
## Database schema (dot)

```dot
digraph G{
vcffile;
sample;
sample2file;
allele;
filter;
chromosome;
variant;
variant2alt;
variant2filter;
vepPrediction;
vepPrediction2so;
genotype;
sample2file -> vcffile[label=vcffile_id];
sample2file -> sample[label=sample_id];
filter -> vcffile[label=vcffile_id];
chromosome -> vcffile[label=vcffile_id];
variant -> vcffile[label=vcffile_id];
variant -> chromosome[label=chromosome_id];
variant -> allele[label=ref_id];
variant2alt -> variant[label=variant_id];
variant2alt -> allele[label=alt_id];
variant2filter -> variant[label=variant_id];
variant2filter -> filter[label=filter_id];
vepPrediction -> variant[label=variant_id];
vepPrediction2so -> vepPrediction[label=vepPrediction_id];
genotype -> variant[label=variant_id];
genotype -> sample[label=sample_id];
genotype -> allele[label=a1_id];
genotype -> allele[label=a2_id];
}

```
END_DOC
*/
@Program(name="vcf2sql",
		description="Generate the SQL code to insert a VCF into mysql",
		keywords={"vcf","sql"}
		)
public class VcfToSql extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfToSql.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-s","--schema"},description="Print Schema")
	private boolean print_schema = false;
	@Parameter(names={"-d","--drop"},description="Add Drop Tables Statement")
	private boolean drop_tables = false;
	@Parameter(names={"-n","--noinfo"},description="ignore INFO column")
	private boolean ignore_info = false;
	@Parameter(names={"-f","--nofilter"},description="ignore FILTER column")
	private boolean ignore_filter = false;
    private PrintWriter outputWriter =null;
    
    private class SelectStmt
    	{
    	String sql;
    	
    	SelectStmt(final Table t,final String field,final Object o)
			{
    		final Column c=t.getColumnByBame(field);
			this.sql="SELECT id from "+t.getAntiquote()+" where "+c.getAntiquote()+"="+c.escape(o);
			}
    	SelectStmt(final Table t)
			{
			this.sql="SELECT max(id) from "+t.getAntiquote();
			}
    	
    	@Override
    	public String toString() {
    		return sql;
    		}
    	}
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
    	public String getAntiquote()
    		{
			return "`"+getName()+"`";
			}
    	
    	}
    
    private class ColumnBuilder
    	{
    	private Table _foreignTable=null;
    	private boolean _pkey=false;
    	private boolean _indexed=false;
    	private boolean _unique=false;
    	private int _length=0;
    	
    	private String _name=null;
    	private Class<?> _class=String.class;
    	private boolean _nilleable=false;
    	ColumnBuilder primaryKey() { this._pkey=true; return this;}
    	ColumnBuilder name(String v) { this._name=v; return this;}
    	ColumnBuilder type(Class<?> v) { this._class=v; return this;}
    	ColumnBuilder foreignKey(Table v) { this._foreignTable=v;return this;}
    	ColumnBuilder foreignKey(Table v,String name) {return foreignKey(v).name(name);}
    	ColumnBuilder nilleable() { this._nilleable=true; return this;}
    	ColumnBuilder indexed() { this._indexed=true; return this;}
    	ColumnBuilder uniq() { this._unique=true; return this;}
    	ColumnBuilder length(int v) { this._length=v; return this;}
    	
    	Column make()
    			{
    			Column c=null;
    			if(_foreignTable!=null)
    				{
    				c = new ForeignKey(this._foreignTable);
    				if(_name!=null) c.name=_name;
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
    			c.unique=this._unique;
    			c.indexed=this._indexed;
    			c.nilleable=this._nilleable;
    			c.maxLength=this._length;
    			if(_pkey) c.nilleable=true;
    			return c;
    			}
    	}
    
    
    private class TableBuilder
    	{
    	boolean _insertIgnore=false;
    	String _name=null;
    	List<Column> _columns=new ArrayList<>();
    	TableBuilder name(String v) { this._name=v; return this;}
    	TableBuilder insertIgnore() { this._insertIgnore=true; return this;}
    	TableBuilder columns(Column...cols) { this._columns.addAll(Arrays.asList(cols)); return this;}
    	
    	
    	
    	public Table make()
    		{
    		Table t=new Table(_name,this._columns);
    		t.insertIgnore=_insertIgnore;
    		return t;
    		}
    	}
    
    private abstract class Column
    	extends AbstractComponent
    	{
    	boolean indexed=false;
    	boolean unique=false;
    	boolean nilleable=false;
    	int maxLength=0;
    	Table table;
    	Column(String name)
			{
			super(name);
			}
    	Object escape(final Object o)
    		{
    		if(o==null)
    			{
    			if(!this.nilleable) throw new RuntimeException("column "+ this.table.getName()+"."+this.getName()+" : not set as nilleable");
    			return "NULL";
    			}
    		else
    			{
    			return o;
    			}
    		}

    	public  void createIndex(PrintWriter pw)
    		{
    		if(unique)
    			{
    			pw.print(",CONSTRAINT "+table.getName()+"_"+getName()+"_uniq UNIQUE("+getAntiquote()+")");
    			}
    		else if( indexed)
    			{
    			pw.print(",INDEX("+getAntiquote()+")");
    			}
    		}
    	public abstract void createColumn(PrintWriter pw);
    	
    	
    	}
    
    private class LongColumn
	extends Column
		{
		LongColumn(String name)
			{
			super(name);
			}
		@Override
		public void createColumn(PrintWriter pw)
			{
			pw.print(getAntiquote()+" INT "+(nilleable?"":" NOT ")+"NULL");
			createIndex(pw);
			}
		}

    
    private class ForeignKey
	extends LongColumn
		{
		Table referencesTable;
    	public ForeignKey(Table t)
    		{
			super(t.getName()+"_id");
			this.referencesTable=t;
    		}
    	
    	
    	@Override
    	public void createIndex(PrintWriter pw) {
    		pw.print( ",FOREIGN KEY ("+getAntiquote()+")  " +
    				"REFERENCES "+referencesTable.getAntiquote()+"(id)"
    				); 
    		}
    	
    	@Override
    	public void createColumn(PrintWriter pw)
    		{
    		pw.print( getAntiquote()+" INT "+(nilleable?"":"NOT")+" NULL");
    		createIndex(pw);
    		}


		}
    
    private class PrimaryKey
   	extends LongColumn
   		{   		
       	public PrimaryKey()
       		{
   			super("id");
       		}
       	       	
       	@Override
       	public void createIndex(PrintWriter pw) {
       		pw.print(",PRIMARY KEY ("+getAntiquote()+")"  ); 
       		}
       	
       	@Override
    	public void createColumn(PrintWriter pw)
    		{
       		pw.print( getAntiquote()+" INT NOT NULL AUTO_INCREMENT");
       		createIndex(pw);
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
		public void createColumn(PrintWriter pw)
			{
			pw.print(  getAntiquote()+" INT "+(nilleable?"":" NOT ")+"NULL" );
			createIndex(pw);
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
		public void createColumn(PrintWriter pw)
			{
			pw.print(  getAntiquote()+" DOUBLE "+(nilleable?"":" NOT ")+"NULL");
			createIndex(pw);
			}

		}

    
    private class StringColumn
    	extends Column
    	{
    	StringColumn(String name)
			{
			super(name);
			}
    	@Override
    	Object escape(Object o)
			{
    		if(o==null)
				{
				return super.escape(o);
				}
			else
				{
				String s=String.valueOf(o);
				if(s.length() >this.maxLength)
					{
					throw new RuntimeException("string length("+s+") greater  than "+this.maxLength+" L="+s.length()+" . Update source code for "+getAntiquote()+" "+table.getName());
					}
				StringBuilder sb=new StringBuilder(s.length()+2);
				sb.append("\"");
				for(int i=0;i< s.length();++i)
					{
					switch(s.charAt(i))
						{	
						case '\"': sb.append("\\\""); break;
						default: sb.append(s.charAt(i)); break;
						}
					}
				sb.append("\"");
				return sb.toString();
				}
			}

    	
		@Override
		public void createColumn(PrintWriter pw)
			{
			pw.print(  getAntiquote()+" VARCHAR("+(maxLength+1)+") "+
					(nilleable?"":" NOT ")+"NULL" );
			createIndex(pw);
			}

    	}
    
    private class Table
    	extends AbstractComponent
    	{
    	boolean insertIgnore=false;
    	List<Column> columns=new ArrayList<>();
    	Map<String,Column> name2column=new HashMap<String,Column>();
    	
    	Table(String name,List<Column> columns)
			{
			super(name);
			this.columns.addAll(columns);
			for(Column c:columns)
				{
				c.table=this;
				this.name2column.put(c.getName(), c);
				}
			}
    	
    	public Column getColumnByBame(String s)
    		{
    		for(Column c:this.columns) if(s.equals(c.getName())) return c;
    		throw new RuntimeException("Cannot find col \""+s +"\" in "+getName());
    		}
    	
    	void insert(PrintWriter pw,Object...row)
    		{
			pw.print("INSERT "+(insertIgnore?"IGNORE":"")+" INTO ");
			pw.print(getAntiquote());
			pw.print("(");
			

    		for(int i=0;i < this.columns.size();++i)
				{
				if(i>0) pw.print(',');
				pw.print(this.columns.get(i).getAntiquote());
				}
    		pw.print(") VALUES (");
    		
			for(int i=0;i < this.columns.size();++i)
    			{
				Column c= this.columns.get(i);
    			if(i>0 ) pw.print(',');
    			if(row[i] instanceof SelectStmt)
    				{
    				pw.print("(");
    				pw.print(SelectStmt.class.cast(row[i]));
    				pw.print(")");
    				}
    			else
    				{
    				pw.print(c.escape(row[i]));
    				}
    			}
			pw.println(");");
    			
    		
    		}
    	
    	public void createTable(PrintWriter pw)
    		{
    		pw.println("CREATE TABLE IF NOT EXISTS "+getAntiquote()+" (" );
    		for(int i=0;i< columns.size();++i)
    			{
    			if(i>0) pw.println(",");
    			columns.get(i).createColumn(pw);
    			}
    		pw.println("\n) ENGINE=InnoDB, DEFAULT CHARSET=utf8 ;");
    		}
    	
    	
    	
    	
    	}
    
    private int MAX_ALLELE_LENGTH=250;
    
    private Table vcfFileTable = new TableBuilder().name("vcffile").columns(
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().name("file").length(250).uniq().make()
    		).make();
    
    private Table sampleTable = new TableBuilder().name("sample").columns(
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().name("name").length(50).uniq().make()
    		).insertIgnore().make();
    
    private Table sample2fileTable = new TableBuilder().name("sample2file").columns(
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(vcfFileTable).make(),
    		new ColumnBuilder().foreignKey(sampleTable).make()
    		).make();

    
    private Table filterTable = new TableBuilder().name("filter").columns(
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(vcfFileTable).make(),
    		new ColumnBuilder().name("name").length(50).make(),
    		new ColumnBuilder().name("description").length(250).make()
    		).make();

    private Table chromosomeTable = new TableBuilder().name("chromosome").columns(
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(vcfFileTable).make(),
    		new ColumnBuilder().name("name").length(50).uniq().make(),
    		new ColumnBuilder().name("chromLength").type(Integer.class).make()
    		).make();
    
    private Table alleleTable = new TableBuilder().name("allele").columns(
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().name("bases").length(MAX_ALLELE_LENGTH).uniq().make()
    		).insertIgnore().make(); 

    
    private Table variantTable = new TableBuilder().name("variant").columns(
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(vcfFileTable).make(),
    		new ColumnBuilder().name("index_in_file").type(Integer.class).make(),
    		new ColumnBuilder().foreignKey(chromosomeTable).make(),
    		new ColumnBuilder().name("pos").indexed().type(Integer.class).make(),
    		new ColumnBuilder().name("rsid").nilleable().length(15).indexed().make(),
    		new ColumnBuilder().foreignKey(alleleTable,"ref_id").make(),
    		new ColumnBuilder().name("qual").nilleable().type(Double.class).make()
    		).make(); 
    
    private Table variant2altTable = new TableBuilder().name("variant2alt").columns(
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(variantTable).make(),
    		new ColumnBuilder().foreignKey(alleleTable,"alt_id").make()
    		).make(); 
    
    private Table variant2filters = new TableBuilder().name("variant2filter").columns(
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(variantTable).make(),
    		new ColumnBuilder().foreignKey(filterTable).make()
    		).make(); 

    private Table vepPrediction = new TableBuilder().name("vepPrediction").columns(
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(variantTable).make(),
    		new ColumnBuilder().name("ensGene").nilleable().length(20).make(),
    		new ColumnBuilder().name("ensTranscript").nilleable().length(20).make(),
    		new ColumnBuilder().name("ensProtein").nilleable().length(20).make(),
    		new ColumnBuilder().name("geneSymbol").nilleable().length(20).make()
    		).make(); 
    
    private Table soTermTable = new TableBuilder().name("soTerm").columns(
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().name("acn").uniq().length(11).make(),
    		new ColumnBuilder().name("description").length(255).make()
    		).insertIgnore().make(); 

    
    private Table vepPrediction2so = new TableBuilder().name("vepPrediction2so").columns(
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(vepPrediction).make(),
    		new ColumnBuilder().foreignKey(soTermTable).make()
    		).make(); 

    
    private Table genotypeTable = new TableBuilder().name("genotype").columns(
    		new ColumnBuilder().primaryKey().make(),
    		new ColumnBuilder().foreignKey(variantTable).make(),
    		new ColumnBuilder().foreignKey(sampleTable).make(),
    		new ColumnBuilder().foreignKey(alleleTable,"a1_id").nilleable().make(),
    		new ColumnBuilder().foreignKey(alleleTable,"a2_id").nilleable().make(),
    		new ColumnBuilder().name("dp").nilleable().type(Integer.class).make(),
    		new ColumnBuilder().name("gq").nilleable().type(Double.class).make()
    		).make(); 

    private Table all_tables[]=new Table[]{
    	   vcfFileTable,sampleTable,sample2fileTable,
    	   alleleTable,filterTable,chromosomeTable,
    	   variantTable,variant2altTable,variant2filters,
    	   vepPrediction,
    	   soTermTable,
    	   vepPrediction2so,
    	   genotypeTable
    	};
    
    
    public VcfToSql()
    	{
    	
    	}
    
	
	private void read(File filename)
		throws IOException
		{

		/* insert ATGC */
		this.alleleTable.insert(outputWriter,null,"A");
		this.alleleTable.insert(outputWriter,null,"C");
		this.alleleTable.insert(outputWriter,null,"G");
		this.alleleTable.insert(outputWriter,null,"T");

		
		/* insert this sample */
		this.vcfFileTable.insert(outputWriter,null,filename);
		final SelectStmt vcffile_id = new SelectStmt(this.vcfFileTable);
		
		final Map<String,SelectStmt> sample2sampleid = new HashMap<String,SelectStmt>();
		final Map<String,SelectStmt> filter2filterid = new HashMap<String,SelectStmt>();
		final Map<String,SelectStmt> chrom2chromId = new HashMap<String,SelectStmt>();
		
		final VCFIterator r=VCFUtils.createVCFIteratorFromFile(filename);
		final VCFHeader header=r.getHeader();
		
		/* parse samples */
		for(final String sampleName:header.getSampleNamesInOrder())
			{
			this.sampleTable.insert(outputWriter,null,sampleName);
			SelectStmt sample_id = new SelectStmt(this.sampleTable, "name", sampleName);
			sample2sampleid.put(sampleName,sample_id);
			
			this.sample2fileTable.insert(outputWriter,null,vcffile_id,sample_id);
			}
		
		/* parse filters */
		for(final VCFFilterHeaderLine filter:header.getFilterLines())
			{
			this.filterTable.insert(
					outputWriter,
					null,
					vcffile_id,
					filter.getID(),
					filter.getValue()
					);
			filter2filterid.put(filter.getID(), new SelectStmt(this.filterTable, "name", filter.getID()));
			}
		filter2filterid.put(VCFConstants.PASSES_FILTERS_v4, new SelectStmt(this.filterTable, "name", VCFConstants.PASSES_FILTERS_v4));

		
		final SAMSequenceDictionary dict= header.getSequenceDictionary();
		if(dict==null)
			{
			throw new RuntimeException("dictionary missing in VCF");
			}
		/* parse sequence dict */
		for(final SAMSequenceRecord ssr: dict.getSequences())
			{
			this.chromosomeTable.insert(
					outputWriter,
					null,
					vcffile_id,
					ssr.getSequenceName(),
					ssr.getSequenceLength()
					);
			chrom2chromId.put(ssr.getSequenceName(), new SelectStmt(this.chromosomeTable,"name",ssr.getSequenceName()));
			}
		
		VepPredictionParser vepPredictionParser=new VepPredictionParserFactory(header).get();
		
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
		int nVariants=0;
		while(r.hasNext())
			{
			if(this.outputWriter.checkError()) break;

			
			VariantContext var= progress.watch(r.next());
			++nVariants;
			/* insert ref allele */
			this.alleleTable.insert(outputWriter, null,var.getReference().getBaseString());
			
			/* insert variant */
			 this.variantTable.insert(
				outputWriter,
				null,
				vcffile_id,
				nVariants,
				chrom2chromId.get(var.getContig()),
				var.getStart(),
				(var.hasID()?var.getID():null),
				new SelectStmt(this.alleleTable, "bases", var.getReference().getBaseString()),
				(var.hasLog10PError()?var.getPhredScaledQual():null)
				);
			
			SelectStmt variant_id = new SelectStmt(variantTable);
			 
			 
			/* insert alternate alleles */
			for(Allele alt: var.getAlternateAlleles())
				{
				/* insert alt allele */
				this.alleleTable.insert(outputWriter, null,alt.getBaseString());

				
				this.variant2altTable.insert(
					outputWriter,
					null,
					variant_id,
					new SelectStmt(this.alleleTable, "bases", alt.getBaseString())
					);
				}

			/* insert filters */
			for(final String filter:var.getFilters())
				{
				if(filter2filterid.get(filter)==null)
					{
					throw new IOException("VCF Error: filter "+filter+" is not defined in the VCF header.");
					}
				this.variant2filters.insert(
					outputWriter,
					null,
					variant_id,
					filter2filterid.get(filter)
					);
				}
			
			if(!this.ignore_info)
				{
				for(final VepPrediction pred: vepPredictionParser.getPredictions(var))
					{
					/*
					vepPrediction.insert(
							outputWriter,
							null,
							variant_id,
							pred.getEnsemblGene(),
							pred.getEnsemblTranscript(),
							pred.getEnsemblProtein(),
							pred.getSymbol()
							);
					SelectStmt pred_id = new SelectStmt(vepPrediction);
			
					for(SequenceOntologyTree.Term t: pred.getSOTerms())
						{
						String term=t.getAcn().replace(':', '_');
						soTermTable.insert(
								outputWriter,
								null,
								term,
								t.getAcn()
								);//for bioportal compatibility
						SelectStmt term_id = new SelectStmt(soTermTable,"acn",term);
						
						vepPrediction2so.insert(
							outputWriter,
							null,
							pred_id,
							term_id
							);
						}
					*/
					}
				}
			
			/* insert genotypes */
			for(final String sampleName: sample2sampleid.keySet())
				{
				final Genotype g= var.getGenotype(sampleName);
				
				if(!g.isAvailable() || g.isNoCall()) continue;
				genotypeTable.insert(
						outputWriter,
						null,
						variant_id,
						sample2sampleid.get(sampleName),
						g.isCalled()?new SelectStmt(this.alleleTable, "bases", g.getAllele(0).getBaseString()):null,
						g.isCalled()?new SelectStmt(this.alleleTable, "bases", g.getAllele(1).getBaseString()):null,
						g.hasDP()?g.getDP():null,
						g.hasGQ()?g.getGQ():null	
						);
				}
			
			}
		r.close();
		}
	
    
	@Override
	public int doWork(final List<String> args) {
		try
			{
			if(this.print_schema)
				{
				this.outputWriter =  this.openFileOrStdoutAsPrintWriter(this.outputFile);
				
				this.outputWriter.println("digraph G{");
				for(int i=0;i< this.all_tables.length;++i)
					{
					this.outputWriter.println(this.all_tables[i].getName()+";");
					}
				for(int i=0;i< this.all_tables.length;++i)
					{
					for(Column c2:this.all_tables[i].columns)
						{
						if(!(c2 instanceof ForeignKey)) continue;
						ForeignKey fk=ForeignKey.class.cast(c2);
						this.outputWriter.println(
								this.all_tables[i].getName()+
								" -> "+
								fk.referencesTable.getName()+
								"[label="+fk.name+"];"
								);
						}
					}
				this.outputWriter.println("}");
				this.outputWriter.flush();
				this.outputWriter.close();
				return RETURN_OK;
				}
			
			//final String inputName=;
			final File filename=new File( oneAndOnlyOneFile(args));
			
			this.outputWriter =  this.openFileOrStdoutAsPrintWriter(this.outputFile);
			
			if(this.drop_tables)
				{
				for(int i=this.all_tables.length-1;i>=0;--i)
					{
		    		this.outputWriter.println("DROP TABLE IF EXISTS "+all_tables[i].getAntiquote()+";");
	
					}
				}
			this.outputWriter.println("START TRANSACTION;");
			this.outputWriter.println("SET autocommit=0;");
			for(Table t:this.all_tables)
				{
				t.createTable(outputWriter);
				}
			
			read(filename);
			
			this.outputWriter.println("COMMIT;");
			this.outputWriter.flush();
			this.outputWriter.close();
			this.outputWriter=null;
			LOG.info("done");
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.outputWriter);
			}
		}
    


	public static void main(final String[] args)
		{
		new VcfToSql().instanceMainWithExit(args);
		}
	}
