/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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

*/
package com.github.lindenb.jvarkit.tools.vcfucsc;

import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter.OnNotFound;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
/**
BEGIN_DOC

## History

20180206: faster creating a prepared statement for each bin.size. fix chromContig

## Example


```
java -jar dist/vcfucsc.jar --table snp142 -e '${name}' input.vcf
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr3	124290753	.	G	C	579.77	.	UCSC_HG19_SNP142=rs145115089
chr3	124290943	.	A	G	491.77	.	UCSC_HG19_SNP142=rs7372055
chr3	124291069	.	G	A	266.77	.	UCSC_HG19_SNP142=rs7373767
chr3	124291171	.	C	CA	240.73	.	.
chr3	124291245	.	A	G	563.77	.	UCSC_HG19_SNP142=rs12695439
chr3	124291351	.	A	G	194.77	.	UCSC_HG19_SNP142=rs7613600
chr3	124291416	.	G	T	308.77	.	UCSC_HG19_SNP142=rs73189597
chr3	124291579	.	T	C	375.77	.	UCSC_HG19_SNP142=rs7649882
```
## Example

```
 java -jar dist/vcfucsc.jar --table vistaEnhancers  --tag VISTAENHANCERS -x 1000 -e '${chromStart}|${chromEnd}|${name}|${score}' input.vcf

```

END_DOC
 */

@Program(
		name="vcfucsc",
		description="annotate an VCF with mysql UCSC data",
		keywords={"ucsc","mysql","vcf"}
		)
public class VcfUcsc extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfUcsc.class).make();

	private static abstract class Printf
		{
		Printf next=null;
		abstract void eval(final ResultSet row,final StringBuilder sb) throws SQLException;
		void append(final Printf n) {if( this.next==null) {this.next=n;} else this.next.append(n);}
		boolean hasColumnIndex(){ return (this.next==null?false:this.next.hasColumnIndex());}
		}
	private static class PrintfText extends Printf
		{
		final String content;
		PrintfText(final String text) {this.content=text;}
		@Override void eval(final ResultSet row,final StringBuilder sb) throws SQLException
			{
			sb.append(content);
			if(next!=null) next.eval(row, sb);
			}
		}

	private static class PrintfColumnIndex extends Printf
		{
		final int column1;
		PrintfColumnIndex(int c) {this.column1=c;}
		@Override  void eval(final ResultSet row,final StringBuilder sb) throws SQLException
			{
			if(column1>0 || column1<=row.getMetaData().getColumnCount())
				{
				sb.append(row.getString(column1));
				}
			
			if(next!=null) next.eval(row, sb);
			}
		boolean hasColumnIndex() {
			return true;
			}
		}

	private static class PrintfColumnLabel extends Printf
		{
		final String label;
		PrintfColumnLabel(String text) {this.label=text;}
		@Override  void eval(ResultSet row,final StringBuilder sb) throws SQLException
			{
			if(!label.isEmpty())
				{
				sb.append(row.getString(label));
				}
			
			if(next!=null) next.eval(row, sb);
			}
		}

	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-e","--expression"},description="expression string.",required=true)
	private String expressionStr="";
	private Printf expression=new PrintfText("");
	
	@Parameter(names={"-D","--database"},description="database name")
	private String database="hg19";
	@Parameter(names={"-T","-t","--table"},description="table name",required=true)
	private String table=null;
	@Parameter(names={"-tag","--tag"},description="INFO tag.")
	private String infoTag=null;
	@Parameter(names={"-jdbc","--jdbc"},description="Java Database Connectivity (JDBC) URI")
	private String jdbcuri="jdbc:mysql://genome-mysql.cse.ucsc.edu";
	@Parameter(names={"-x","--extend"},description="Extend variant coordinates by 'x' bases.")
	private int extend_bases = 0;

	
	private Connection connection=null;
	private boolean has_bin_column=false;
	private String chromColumn=null;
	private String startColumn=null;
	private String endColumn=null;
	private final ContigNameConverter contigNameConverter = ContigNameConverter.createConvertToUcsc().setOnNotFound(OnNotFound.RETURN_ORIGINAL);
	private final Set<String> userColumnNames = new LinkedHashSet<>();
	 
	private void select(final Set<String> atts,final PreparedStatement pstmt) throws SQLException
		{
		final ResultSet row=pstmt.executeQuery();
		while(row.next())
			{
			final StringBuilder sb=new StringBuilder();
			this.expression.eval(row,sb);
			final  String s=sb.toString();
			if(StringUtil.isBlank(s)) continue;
			atts.add(s);
			}
		row.close();
		}
	
	 private static List<Integer> reg2bins(final int beg, final int _end) {
	        int k, end = _end;
	        if (beg >= end) return Collections.emptyList();
	        if (end >= 1 << 29) end = 1 << 29;
	        --end;
	        final List<Integer> list = new ArrayList<>();
	        list.add(0);
	        for (k = 1 + (beg >> 26); k <= 1 + (end >> 26); ++k) list.add(k);
	        for (k = 9 + (beg >> 23); k <= 9 + (end >> 23); ++k) list.add(k);
	        for (k = 73 + (beg >> 20); k <= 73 + (end >> 20); ++k) list.add(k);
	        for (k = 585 + (beg >> 17); k <= 585 + (end >> 17); ++k) list.add(k);
	        for (k = 4681 + (beg >> 14); k <= 4681 + (end >> 14); ++k) list.add(k);
	        return list;
	    }
	 
	 private void initPstmt(final PreparedStatement pstmt,final String contig,int start0,int end0) throws java.sql.SQLException
	 	{
		pstmt.setString(1, this.contigNameConverter.apply(contig));
		pstmt.setInt(2, start0 ) ;
		pstmt.setInt(3, end0);
	 	}
	 
	 private PreparedStatement createPreparedStatement(int nBins) throws SQLException
	 	{
		final StringBuilder b=new StringBuilder("select ");
		if(this.expression.hasColumnIndex())
			{
			b.append(" * ");
			}
		else
			{
			b.append(String.join(",",this.userColumnNames));
			}
		
		b.append(" from "+database+"."+table+" where ");
		b.append(chromColumn).append("=? and NOT(");//contig
		b.append(endColumn).append("<=? or ");//start0
		b.append(startColumn).append(">=? ) ");//end0
		
		if(nBins>0)
			{
			b.append(" and (");
			for(int i=0;i< nBins;i++)
				{
				if(i>0) b.append(" OR ");
				b.append("bin=?");
				}
			b.append(" )");
			}
		
		return this.connection.prepareStatement(b.toString());
		}
	
	@Override
	protected int doVcfToVcf(
			final String inputName, 
			final VcfIterator in, 
			final VariantContextWriter out
			) 
		{
		final String TAG;
		
		if(StringUtil.isBlank(this.infoTag))
			{
			TAG = "UCSC_"+database.toUpperCase()+"_"+table.toUpperCase();
			}
		else
			{
			TAG= this.infoTag;
			}
		ResultSet row=null;
		VCFHeader header=in.getHeader();
		
		final VCFHeader h2=new VCFHeader(header);
		h2.addMetaDataLine(new VCFInfoHeaderLine(
				TAG,
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,
				database+"."+table)
				);
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		out.writeHeader(h2);
		
		final Map<Integer, PreparedStatement> bin2pstmt = new HashMap<>();
		try	{			
			if(!this.has_bin_column)
				{
				bin2pstmt.put(0, createPreparedStatement(0));
				}
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header).logger(LOG);
			while(in.hasNext())
				{
				final VariantContext ctx= progress.watch(in.next());
				int start0 ,end0;
				
				if(ctx.isIndel()) //mutation starts *after* the base
					{
					start0 = ctx.getStart();
					end0 = ctx.getEnd();
					}
				else
					{
					start0 = ctx.getStart() -1;
					end0 = ctx.getEnd();
					}
				// extends left/right
				start0 = Math.max(0, start0-this.extend_bases);
				end0 += this.extend_bases;
				
				final Set<String> atts=new HashSet<String>();
				if(this.has_bin_column)
					{
					final List<Integer> binList = reg2bins(start0, end0);
					PreparedStatement pstmt = bin2pstmt.get(binList.size());
					if(pstmt==null) {
						LOG.debug("create prepared statemement for bin.size="+binList.size()+"["+start0+":"+end0+"]");
						pstmt = createPreparedStatement(binList.size());
						bin2pstmt.put(binList.size(), pstmt);
						}
					initPstmt(pstmt,ctx.getContig(),start0,end0);
					for(int x=0;x< binList.size();++x)
						{
						pstmt.setInt(4+x, binList.get(x));
						}
					
					select(atts,pstmt);
					}
				else
					{
					final PreparedStatement pstmt = bin2pstmt.get(0);//alread defined
					initPstmt(pstmt,ctx.getContig(),start0,end0);
					select(atts,pstmt);
					}
				if(atts.isEmpty())
					{
					out.add(ctx);
					continue;
					}
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.attribute(TAG,atts.toArray());
				out.add(vcb.make());
				}
			progress.finish();
			return 0;
			}
		catch(final SQLException err)
			{
			LOG.error(err);
			throw new RuntimeException("SQLError", err);
			}
		finally
			{
			for(PreparedStatement pstmt: bin2pstmt.values()) 
				{
				CloserUtil.close(pstmt);
				}
			bin2pstmt.clear();
			CloserUtil.close(row);
			}
		}

	@Override
	public int doWork(final List<String> args) {
		int max_column_index=0;
		try
			{
			
			if(StringUtil.isBlank(this.expressionStr))
				{
				LOG.error("expression missing");
				return -1;
				}
			if(StringUtil.isBlank(this.table))
				{
				LOG.error("Table undefined.");
				return -1;
				}

			
			String s = this.expressionStr;
			while(!s.isEmpty())
				{
				int b=s.indexOf("${");
				if(b==-1)
					{
					this.expression.append(new PrintfText(s));
					break;
					}
				int e=s.indexOf("}",b);
				if(e==-1) 
					{
					LOG.error("Cannot find end of expression in "+this.expressionStr);
					return -1;
					}
				if(b>0) this.expression.append(new PrintfText(s.substring(0,b)));
				final String between = s.substring(b+2,e).trim();
				if(between.isEmpty()) 
					{
					LOG.error("bad expression in "+this.expressionStr);
					return -1;
					}
				if(between.matches("[0-9]+"))
					{
					int colIndex = Integer.parseInt(between);
					max_column_index = Math.max(max_column_index, colIndex);
					this.expression.append(new PrintfColumnIndex(colIndex));
					}
				else
					{
					this.userColumnNames.add(between);
					this.expression.append(new PrintfColumnLabel(between));
					}
				s=s.substring(e+1);
				}
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
					
			
		
		if(this.expression==null)
			{
			LOG.error("Expression undefined.");
			return -1;
			}
		try
			{
			LOG.info("Getting jdbc-driver");
			Class.forName("com.mysql.jdbc.Driver");
			
			this.connection=DriverManager.getConnection(
					jdbcuri+"/"+database+"?user=genome&password=");
			LOG.info("Getting jdbc-driver: Done.");
			final Statement stmt=this.connection.createStatement();
			final ResultSet row=stmt.executeQuery("describe "+this.database+"."+this.table);
			final Set<String> cols=new HashSet<String>();
			while(row.next())
				{
				final String colName=row.getString("Field");
				if(StringUtil.isBlank(colName)) {
					LOG.error("empty field in "+this.database+"."+this.table);
					return -1;
					}
				cols.add(colName);
				}
			row.close();
			stmt.close();
			this.has_bin_column=cols.contains("bin");
			
			if(max_column_index> cols.size())
				{
				LOG.error("No column index["+max_column_index+"] for "+cols+" N="+cols.size());
				return -1;
				}
			
			for(final String userColumnName: this.userColumnNames)
				{
				if(!cols.contains(userColumnName)) {
					LOG.error("No column '"+userColumnName+"' in "+cols);
					return -1;
					}
				}
			
			for(final String col:new String[]{"chrom"})
				{
				if(this.chromColumn==null && cols.contains(col))
					{
					this.chromColumn=col;
					}
				}
			if(this.chromColumn==null)
				{
				LOG.error("cannot find 'chrom' in the columns of '"+this.database+"."+this.table+ "' : "+cols);
				return -1;
				}
			
			for(final String col:new String[]{"txStart","cdsStart","chromStart"})
				{
				if(this.startColumn==null && cols.contains(col))
					{
					this.startColumn=col;
					}
				}
			if(this.startColumn==null)
				{
				LOG.error("cannot find startColumn in "+cols);
				return -1;
				}
			for(final String col:new String[]{"txEnd","cdsEnd","chromEnd"})
				{
				if(this.endColumn==null && cols.contains(col))
					{
					this.endColumn=col;
					}
				}
			if(this.endColumn==null)
				{
				LOG.error("cannot find endColumn in "+cols);
				return -1;
				}
			return doVcfToVcf(args,this.outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.connection);
			}
		}

	public static void main(String[] args)
		{
		new VcfUcsc().instanceMainWithExit(args);
		}

	}
