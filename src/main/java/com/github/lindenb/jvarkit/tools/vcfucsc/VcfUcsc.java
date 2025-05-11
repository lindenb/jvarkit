/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.function.Predicate;

import org.apache.commons.jexl2.JexlContext;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jexl.BaseJexlHandler;
import com.github.lindenb.jvarkit.jexl.JexlPredicate;
import com.github.lindenb.jvarkit.jexl.JexlToString;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

@Program(
		name="vcfucsc",
		description="annotate an VCF with mysql UCSC data",
		keywords={"ucsc","mysql","vcf"},
		creationDate="20160105",
		modificationDate="20210201"
		)
public class VcfUcsc extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.of(VcfUcsc.class);

	
	private static class ResultSetJexlContext implements JexlContext {
		final VariantContext variant;
		final ResultSet row;
		final ResultSetMetaData metadata;
		ResultSetJexlContext(final VariantContext variant,final ResultSet row,final ResultSetMetaData metadata) {
			this.variant  = variant;
			this.row  = row;
			this.metadata = metadata;
			}
		
		@Override
		public Object get(final String arg) {
			if(arg.equals("variant")) return this.variant;
			if(arg.equals("row")) return this.row;
			if(arg.equals("meta")) return this.metadata;
			try {
				return this.row.getObject(arg);
			} catch (final SQLException e) {
				throw new RuntimeException(e);
				}
			}
		@Override
		public void set(String arg0, Object arg1) {
			throw new UnsupportedOperationException("Cannot set "+arg0);
			}
		@Override
		public boolean has(String arg) {
			if(arg.equals("variant")) return true;
			if(arg.equals("row")) return true;
			if(arg.equals("meta")) return true;
			try {
				for(int i=1;i<=this.metadata.getColumnCount();i++) {
					if(this.metadata.getColumnName(i).equals(arg)) return true;
					}
			} catch (final SQLException e) {
				throw new RuntimeException(e);
				}
			return false;
			}
		}
		
	private Function<JexlContext,String> toStringFunc = ROW->"DEFAULT_OUTPUT";
	private Predicate<JexlContext> acceptRowFunc = ROW->true;

	
	@Parameter(names={"-e","--expression"},description="JEXL expression used to convert a row to String. Empty=default. See the manual. " + BaseJexlHandler.OPT_WHAT_IS_JEXL)
	private String convertToStrExpr ="";
	@Parameter(names={"-a","--accept"},description="JEXL expression used to accept a result set. Must return a boolean. Empty=defaul/accept all. See the manual. " + BaseJexlHandler.OPT_WHAT_IS_JEXL)
	private String acceptExpr="";
	@Parameter(names={"-L","--limit"},description="Limit number of items after filtration. Negative = no-limit")
	private int limit_item_numbers = -1;

	@Parameter(names={"-D","--database"},description="mysql database name.")
	private String database="hg19";
	@Parameter(names={"-T","-t","--table"},description="table name",required=true)
	private String table=null;
	@Parameter(names={"-tag","--tag"},description="INFO tag.")
	private String infoTag=null;
	@Parameter(names={"-jdbc","--jdbc"},description="Java Database Connectivity (JDBC) URI")
	private String jdbcuri="jdbc:mysql://genome-mysql.cse.ucsc.edu";
	@Parameter(names={"-x","--extend"},description="Extend variant coordinates by 'x' bases.")
	private int extend_bases = 0;
	@Parameter(names={"-fi","--filterIn"},description="Set this FILTER if any item is found in the database")
	private String filterIn=null;
	@Parameter(names={"-fo","--filterOut"},description="Set this FILTER if no item is found in the database")
	private String filterOut=null;
	
	private Connection connection=null;
	private boolean has_bin_column=false;
	private String chromColumn=null;
	private String startColumn=null;
	private String endColumn=null;
	private final ContigNameConverter contigNameConverter = ContigNameConverter.createConvertToUcsc();
	 
	private void select(final VariantContext ctx,final Set<String> atts,final PreparedStatement pstmt) throws SQLException
		{
		try( ResultSet row=pstmt.executeQuery()){
			final ResultSetMetaData meta = pstmt.getMetaData();
			while(row.next() && (this.limit_item_numbers <0 || atts.size()< this.limit_item_numbers))
				{
				final ResultSetJexlContext rowCtx = new ResultSetJexlContext(ctx,row,meta);
				if(!this.acceptRowFunc.test(rowCtx)) continue;
				final  String s = this.toStringFunc.apply(rowCtx);
				if(StringUtil.isBlank(s)) continue;
				atts.add(s);
				}
			}
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
		String ctg = this.contigNameConverter.apply(contig);
		if(ctg==null) ctg=contig;
		pstmt.setString(1, ctg);
		pstmt.setInt(2, start0 ) ;
		pstmt.setInt(3, end0);
	 	}
	 
	 private PreparedStatement createPreparedStatement(int nBins) throws SQLException
	 	{
		final StringBuilder b=new StringBuilder("select ");
		b.append(" * ");
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
			final VCFIterator in, 
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
		final VCFHeader header=in.getHeader();
		
		final SAMSequenceDictionary dict = header.getSequenceDictionary();
		if(dict!=null ) {
			if(this.database.equals("hg19") && !SequenceDictionaryUtils.isGRCh37(dict))
				{
				LOG.warn("using hg19 but sequence dictionary doesn't look like hg19 ?");
				}
			else if(this.database.equals("hg38") && !SequenceDictionaryUtils.isGRCh38(dict))
				{
				LOG.warn("using hg38 but sequence dictionary doesn't look like hg38 ?");
				}
			}
		
		
		final VCFHeader h2=new VCFHeader(header);
		h2.addMetaDataLine(new VCFInfoHeaderLine(
				TAG,
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,
				database+"."+table)
				);
		
		if(!StringUtil.isBlank(this.filterIn))
			{
			h2.addMetaDataLine(new VCFFilterHeaderLine(
					this.filterIn,
					"Set by "+this.getClass().getName()
					));
			}
		else if(!StringUtil.isBlank(this.filterOut))
			{
			h2.addMetaDataLine(new VCFFilterHeaderLine(
					this.filterOut,
					"Set by "+this.getClass().getName()
					));
			}
		
		JVarkitVersion.getInstance().addMetaData(this, h2);
		out.writeHeader(h2);
		
		final Map<Integer, PreparedStatement> bin2pstmt = new HashMap<>();
		try	{			
			if(!this.has_bin_column)
				{
				bin2pstmt.put(0, createPreparedStatement(0));
				}
			while(in.hasNext())
				{
				final VariantContext ctx=  in.next();
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
					final PreparedStatement pstmt;
					if(bin2pstmt.containsKey(binList.size())) {
						pstmt =  bin2pstmt.get(binList.size());
						}
					else {
						LOG.debug("create prepared statemement for bin.size="+binList.size()+"["+start0+":"+end0+"]");
						pstmt = createPreparedStatement(binList.size());
						bin2pstmt.put(binList.size(), pstmt);
						}
					initPstmt(pstmt,ctx.getContig(),start0,end0);
					for(int x=0;x< binList.size();++x)
						{
						pstmt.setInt(4+x, binList.get(x));
						}
					
					select(ctx,atts,pstmt);
					}
				else
					{
					final PreparedStatement pstmt = bin2pstmt.get(0);//already defined
					initPstmt(pstmt,ctx.getContig(),start0,end0);
					select(ctx,atts,pstmt);
					}
				if(atts.isEmpty() && StringUtil.isBlank(this.filterIn) && StringUtil.isBlank(this.filterOut))
					{
					out.add(ctx);
					continue;
					}
				
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				if(!StringUtil.isBlank(this.filterIn) &&  !atts.isEmpty())
					{
					vcb.filter(this.filterIn);
					}
				else if(!StringUtil.isBlank(this.filterOut) &&  atts.isEmpty())
					{
					vcb.filter(this.filterOut);
					}
				if(!atts.isEmpty()) {
					vcb.attribute(TAG,new ArrayList<String>(atts));
					}
				out.add(vcb.make());
				}
			return 0;
			}
		catch(final SQLException err)
			{
			LOG.error(err);
			throw new RuntimeException("SQLError", err);
			}
		finally
			{
			for(final PreparedStatement pstmt: bin2pstmt.values()) 
				{
				CloserUtil.close(pstmt);
				}
			bin2pstmt.clear();
			}
		}

	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int beforeVcf() {
		int max_column_index=0;
		try
			{
			
			
			if(StringUtil.isBlank(this.table))
				{
				LOG.error("Table undefined.");
				return -1;
				}

			if(!StringUtil.isBlank(this.filterIn) && !StringUtil.isBlank(this.filterOut))
				{
				LOG.error("both filters in/out defined.");
				return -1;
				}
			
			if(!StringUtil.isBlank(this.acceptExpr)) {
				this.acceptRowFunc = new JexlPredicate(this.acceptExpr);
			}
			if(!StringUtil.isBlank(this.convertToStrExpr)) {
				this.toStringFunc = new JexlToString(this.convertToStrExpr);
			}

			}
		catch(final Throwable err)
			{
			LOG.error(err);
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
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		}
	
	@Override
	protected void afterVcf() {
		try {
			CloserUtil.close(this.connection);
			} catch(final Throwable err) {
			LOG.error(err);
			}
		}

	public static void main(final String[] args)
		{
		new VcfUcsc().instanceMainWithExit(args);
		}

	}
