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
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.util.CloserUtil;

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
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
/**
BEGIN_DOC

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
		abstract void eval(ResultSet row,StringBuilder sb) throws SQLException;
		void append(Printf n) {if( this.next==null) {this.next=n;} else this.next.append(n);}	
		}
	private static class PrintfText extends Printf
		{
		final String content;
		PrintfText(String text) {this.content=text;}
		void eval(ResultSet row,StringBuilder sb) throws SQLException
			{
			sb.append(content);
			if(next!=null) next.eval(row, sb);
			}
		}

	private static class PrintfColumnIndex extends Printf
		{
		final int column1;
		PrintfColumnIndex(int c) {this.column1=c;}
		void eval(ResultSet row,StringBuilder sb) throws SQLException
			{
			if(column1>0 || column1<=row.getMetaData().getColumnCount())
				{
				sb.append(row.getString(column1));
				}
			
			if(next!=null) next.eval(row, sb);
			}

		}

	private static class PrintfColumnLabel extends Printf
		{
		final String label;
		PrintfColumnLabel(String text) {this.label=text;}
		void eval(ResultSet row,StringBuilder sb) throws SQLException
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
	@Parameter(names={"-e","--expression"},description="What should be displayed in the INFO field. expression string using column offsets. or names e.g: \"${chrom}:${2}-${3}\"")
	private String expressionStr="";
	private Printf expression=new PrintfText("");
	
	@Parameter(names={"-D","--database"},description="database name")
	private String database="hg19";
	@Parameter(names={"-T","-t","--table"},description="table name",required=true)
	private String table=null;
	@Parameter(names={"-tag","--tag"},description="tag prefix. Default = UCSC_${db}_${table}")
	private String tag="";
	private Connection connection=null;
	private String jdbcuri="jdbc:mysql://genome-mysql.cse.ucsc.edu";
	private boolean has_bin_column=false;
	private String chromColumn=null;
	private String startColumn=null;
	private String endColumn=null;
	
	private void select(final Set<String> atts,PreparedStatement pstmt) throws SQLException
		{
		final ResultSet row=pstmt.executeQuery();
		while(row.next())
			{
			final StringBuilder sb=new StringBuilder();
			this.expression.eval(row,sb);
			String s=sb.toString();
			if(s.isEmpty()) continue;
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
	
	@Override
	protected int doVcfToVcf(final String inputName, VcfIterator in, VariantContextWriter out) 
		{
		final String TAG= this.tag.isEmpty()?
				"UCSC_"+database.toUpperCase()+"_"+table.toUpperCase():
				this.tag
				;
		PreparedStatement pstmt=null;
		ResultSet row=null;
		final VCFHeader header=in.getHeader();
		if(header.getInfoHeaderLine(TAG)!=null)
			{
			LOG.error("VCFheader already contains INFO="+TAG);
			return -1;
			}
		
		
		final VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFInfoHeaderLine(
				TAG,
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,
				database+"."+table)
				);
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));

		
		out.writeHeader(h2);
		try	{
			final StringBuilder b=new StringBuilder(
					"select * from "+database+"."+table+" where ");
			b.append(chromColumn).append("=? and NOT(");
			b.append(endColumn).append("<=? or ");
			b.append(startColumn).append(">=? ) ");
			if(has_bin_column) b.append(" and bin=?");
			
			pstmt=connection.prepareStatement(b.toString());
			
			while(in.hasNext())
				{
				final VariantContext ctx=in.next();
				final Set<String> atts=new HashSet<String>();
				pstmt.setString(1, ctx.getContig());
				pstmt.setInt(2, ctx.getStart());
				pstmt.setInt(3, ctx.getEnd());
				if(this.has_bin_column)
					{
					for(final Integer biter : reg2bins(ctx.getStart()-1, ctx.getEnd()))
						{
						pstmt.setInt(4, biter);
						select(atts,pstmt);
						}
					}
				else
					{
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
			return 0;
			}
		catch(final SQLException err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(row);
			CloserUtil.close(pstmt);
			}
		}

	@Override
	public int doWork(final List<String> args) {
		
		try
			{
			String s=this.expressionStr;
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
					 this.expression.append(new PrintfColumnIndex(Integer.parseInt(between)));
					}
				else
					{
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
					
			
		
		if(this.table==null)
			{
			LOG.error("Table undefined.");
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
			Statement stmt=this.connection.createStatement();
			ResultSet row=stmt.executeQuery("describe "+database+"."+table);
			Set<String> cols=new HashSet<String>();
			while(row.next())
				{
				String colName=row.getString("Field");
				cols.add(colName);
				}
			row.close();
			stmt.close();
			this.has_bin_column=cols.contains("bin");
			
			
			for(final String col:new String[]{"chrom"})
				{
				if(this.chromColumn==null && cols.contains(col))
					{
					this.chromColumn=col;
					}
				}
			if(chromColumn==null)
				{
				LOG.error("cannot find chromColumn in "+cols);
				return -1;
				}
			
			for(final String col:new String[]{"txStart","cdsStart","chromStart"})
				{
				if(this.startColumn==null && cols.contains(col))
					{
					this.startColumn=col;
					}
				}
			if(startColumn==null)
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
			if(endColumn==null)
				{
				LOG.error("cannot find endColumn in "+cols);
				return -1;
				}
			return doVcfToVcf(args, outputFile);
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
