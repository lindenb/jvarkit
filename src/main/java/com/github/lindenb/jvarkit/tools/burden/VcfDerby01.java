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
* 2015 moving to knime

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.File;
import java.io.FileFilter;
import java.io.PrintWriter;
import java.io.Reader;
import java.sql.Clob;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;
import java.sql.Statement;
import java.sql.Types;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Properties;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;


public class VcfDerby01
	extends AbstractVcfDerby01
	{
	private static int MAX_REF_BASE_LENGTH=50;
	private long ID_GENERATOR = System.currentTimeMillis();
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfDerby01.class);
	private Connection conn=null;
	private static final String VCF_HEADER_FILE_ID="##VcfDerby01VcfId=";
	private static final String VCF_HEADER_FILE_NAME="##VcfDerby01VcfName=";
	public VcfDerby01()
		{
		}
	 
	@Override
	public Collection<Throwable> initializeKnime() {
		if(super.derbyFilePath==null ||super.derbyFilePath.isEmpty())
			{
			return wrapException("Undefined Derby DB Directory: option -"+OPTION_DERBYFILEPATH);
			}
		try {
			Class.forName("org.apache.derby.jdbc.ClientDriver").newInstance();
		} catch(Exception err ){
			LOG.error("Cannot get derby driver");
			return wrapException(err);
		}
		return super.initializeKnime();
	 	}
	
	
	private File getDerbyDirectory() {
		return new File(this.derbyFilePath);
	}
	
	
	private void openDerby() {
		try {
			boolean create;
			final Properties props = new Properties();
			final File derbyDir = getDerbyDirectory();
			LOG.info("open derby :" + getDerbyDirectory());
			if(derbyDir.exists()) {
				if(!derbyDir.isDirectory()) {
					throw new RuntimeIOException("derby database is not a directory : "+derbyDir);
					}
				
				if( derbyDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						if(pathname.isFile()){
							if( pathname.getName().equals("service.properties")) return true;
							if( pathname.getName().equals("README_DO_NOT_TOUCH_FILES.txt")) return true;
							}
						if(pathname.isDirectory()){
							if( pathname.getName().startsWith("log")) return true;
						}
						return false;
					}
				}).length==0) {
					throw new RuntimeIOException("derby database exist but doesn't look like a derby directory : "+derbyDir);
				}
				
				create=false;
				}
			else
				{
				create=true;
				}
			props.setProperty("create", String.valueOf(create));
			this.conn = DriverManager.getConnection("jdbc:derby:"+derbyDir,props);
			
			if(create) {
				final String tableId = "ID INTEGER NOT NULL GENERATED ALWAYS AS IDENTITY (START WITH 1, INCREMENT BY 1) PRIMARY KEY";
				final Statement stmt= this.conn.createStatement();
				final String sqls[]={
						"CREATE TABLE ROWCONTENT("+tableId+",MD5SUM CHAR(32) UNIQUE,CONTENT CLOB,CONTIG VARCHAR(20),FILTERED SMALLINT NOT NULL,START INT,STOP INT,ALLELE_REF VARCHAR("+MAX_REF_BASE_LENGTH+"))",
						"CREATE TABLE VCF("+tableId+",NAME VARCHAR(255))",
						"CREATE TABLE VCFROW("+tableId+",VCF_ID INTEGER CONSTRAINT row2vcf REFERENCES VCF,ROW_ID INTEGER CONSTRAINT row2content REFERENCES ROWCONTENT)"
						};
				for(final String sql:sqls) {
					LOG.warn(sql);
					stmt.execute(sql);
				}
				stmt.close();
			}
			this.conn.setAutoCommit(true);
		} catch (Exception e) {
			CloserUtil.close(this.conn);
			this.conn = null;
			throw new RuntimeException(e);
		}
	}
	
	private static long getLastGeneratedId(final PreparedStatement pstmt) throws SQLException {
			ResultSet keys = null;
			long id = -1L;
			try {
			keys = pstmt.getGeneratedKeys();
		    while (keys.next()) {
		       if(id!=-1L) throw new IllegalStateException();
		       id= keys.getLong(1);
		       }	 
			if(id==-1L) throw new IllegalStateException("No SQL ID was found");
			return id;
			}
			finally {
				CloserUtil.close(keys);
			}
		}
	
	private void closeDerby() {
		CloserUtil.close(this.conn);
		this.conn = null;
		try {
			final Properties props = new Properties();
			props.setProperty("shutdown", "true");
			final File derbyDir = getDerbyDirectory();
			DriverManager.getConnection("jdbc:derby:"+derbyDir,props);
		} catch (Exception e) {
			
			}
		}
	
	private void compress() {
		PreparedStatement pstmt =null;
		
		for(final String sqlStmt :new String[]{
				"call SYSCS_UTIL.SYSCS_COMPRESS_TABLE(?,?, 1)",	
				"call SYSCS_UTIL.SYSCS_INPLACE_COMPRESS_TABLE(?,?, 1, 1, 1)"
		}) {
		try {
		/* compress data */
		pstmt = this.conn.prepareStatement(sqlStmt);
		pstmt.setString(1, "APP");
		for(final String table: new String[]{"ROWCONTENT","VCF","VCFROW"}) {
			LOG.info("compressing "+table+" with "+sqlStmt);
			pstmt.setString(2,table);
			pstmt.execute();
		}
		pstmt.close();
		} catch(Exception err) {
			LOG.warn("Cannot compress", err);
		}
		finally 
		{
			CloserUtil.close(pstmt);
		}
		}

		
		
		
	}
	
	private Collection<Throwable> doCommandDumpAll(){
		Statement stmt = null;
		ResultSet row = null;
		final Set<Long> vcfids = new TreeSet<>();
		if(!getInputFiles().isEmpty()) {
			return wrapException("Too many arguments");
		}
		try {
			stmt = this.conn.createStatement();
			row = stmt.executeQuery("SELECT ID from VCF");
			while(row.next()) {
				vcfids.add(row.getLong(1));
			}
			row.close();row=null;
			stmt.close();stmt=null;
			return dump(vcfids);
		} catch (final Exception e) {
			return wrapException(e);
		} finally {
			CloserUtil.close(row);
			CloserUtil.close(stmt);
		}
	}
	
	private Collection<Throwable> doCommandDumpUniq(){
		if(!getInputFiles().isEmpty()) {
			return wrapException("Too many arguments");
		}
		PreparedStatement pstmt2 = null;
		ResultSet row = null;
		PrintWriter pwOut = null;
		try {
			boolean chrom_line_seen=false;
			for(int side=0;side<2;++side)
				{
				final String sql=(side==0?
						"SELECT ROWCONTENT.CONTENT FROM ROWCONTENT WHERE ROWCONTENT.CONTIG IS NULL ORDER BY ROWCONTENT.ID " :
						"SELECT ROWCONTENT.CONTENT FROM ROWCONTENT WHERE ROWCONTENT.CONTIG IS NOT NULL ORDER BY ROWCONTENT.CONTIG,ROWCONTENT.START,ROWCONTENT.ALLELE_REF "
						);
				LOG.info(sql);
				pstmt2 = this.conn.prepareStatement(sql);
				pwOut = openFileOrStdoutAsPrintWriter();
				row =  pstmt2.executeQuery();
				while(row.next()) {
					final Clob clob = row.getClob(1);
					final Reader r= clob.getCharacterStream();
					if(side==0)
						{
						final String s = IOUtils.copyToString(r);
						if(s.startsWith("#CHROM")) {
							if(chrom_line_seen) {
								r.close();
								continue;
							}
							chrom_line_seen=true;
						} else
							{
							if(chrom_line_seen) {
								r.close();
								continue;
								}
							}
						pwOut.print(s);
						}
					else
						{
						IOUtils.copyTo(r,pwOut);
						}	
					
					
					r.close();
					pwOut.println();
					}
				row.close();
					
				pstmt2.close();pstmt2=null;
				
				pwOut.flush();
				}
			pwOut.close();pwOut=null;
			
			return RETURN_OK;
		} catch (final Exception e) {
			return wrapException(e);
		} finally {
			CloserUtil.close(pwOut);
			CloserUtil.close(row);
			CloserUtil.close(pstmt2);
		}
	}

	private Collection<Throwable> dumpAll(){
		PreparedStatement pstmt2 = null;
		ResultSet row = null;
		PrintWriter pwOut = null;
		int num_vcf_exported=0;
		try {
			pstmt2 = this.conn.prepareStatement(
					"SELECT ROWCONTENT.CONTENT,VCF.ID,VCF.NAME FROM VCF,VCFROW,ROWCONTENT WHERE VCFROW.VCF_ID=VCF.ID AND VCFROW.ROW_ID = ROWCONTENT.ID AND ORDER BY VCF.ID,VCFROW.ID ");
			pwOut = openFileOrStdoutAsPrintWriter();
			row =  pstmt2.executeQuery();
			final String CHROM_prefix="#CHROM\t";
			final StringBuilder chrom_header_line = new StringBuilder(CHROM_prefix.length());
			while(row.next()) {
				final Clob clob = row.getClob(1);
				final Reader r= clob.getCharacterStream();
				/* read the first bytes to check if it's the #CHROM line
				 * if true, add a VCF header line with VCF ID and NAME
				 *  */
				chrom_header_line.setLength(0);
				int c;
				while(chrom_header_line.length()< CHROM_prefix.length() &&
					((c=r.read())!=1) )
					{
					chrom_header_line.append((char)c);
					}
				final String amorce = chrom_header_line.toString();
				if(amorce.equals(CHROM_prefix))
					{
					num_vcf_exported++;
					final long vcf_id = row.getLong(2);
					final String vcfName= row.getString(3);
					pwOut.println(VCF_HEADER_FILE_ID+vcf_id);
					pwOut.println(VCF_HEADER_FILE_NAME+vcfName);
					}
				pwOut.print(amorce);
				IOUtils.copyTo(r,pwOut);
				r.close();
				pwOut.println();
			}
			row.close();
				
			pstmt2.close();
			
			pwOut.flush();
			pwOut.close();
			
			if(num_vcf_exported==0 ) {
				LOG.warn("NO VCF WAS EXPORTED");
				}
			LOG.info("count(VCF) exported "+num_vcf_exported);
			return RETURN_OK;
		} catch (final Exception e) {
			return wrapException(e);
		} finally {
			CloserUtil.close(pwOut);
			CloserUtil.close(row);
			CloserUtil.close(pstmt2);
		}
	}
	

	private Collection<Throwable> dump(final Set<Long> vcfIds){
		final double timeStart = System.currentTimeMillis();

		PreparedStatement pstmt = null;
		PreparedStatement pstmt2 = null;
		ResultSet row = null;
		PrintWriter pwOut = null;
		int num_vcf_exported=0;
		try {
			pstmt = this.conn.prepareStatement("SELECT NAME from VCF where ID=?");
			pstmt2 = this.conn.prepareStatement("SELECT ROWCONTENT.CONTENT FROM VCF,VCFROW,ROWCONTENT WHERE VCFROW.VCF_ID=VCF.ID AND VCFROW.ROW_ID = ROWCONTENT.ID AND VCF.ID=? ORDER BY VCFROW.ID ");
			
			
			pwOut = openFileOrStdoutAsPrintWriter();
			
			
			for(final long vcf_id : vcfIds)
				{
				final double remain = ((((System.currentTimeMillis()-timeStart)/(1+num_vcf_exported)))*(vcfIds.size()-(1+num_vcf_exported)))/1000.0;
				
				LOG.info("Getting VCF "+vcf_id+" "+(num_vcf_exported)+"/"+vcfIds.size() +" . Remains "+(long)remain+" seconds.");
				
				pstmt.setLong(1, vcf_id);
				String vcfName=null;
				row = pstmt.executeQuery();
				while(row.next()) {
					vcfName = row.getString(1);
				}
				row.close();
				if(vcfName==null) 
					{
					return wrapException("Cannot find VCF ID "+vcf_id);
					}
				LOG.info("dumping "+vcfName+" ID:"+vcf_id);
				pstmt2.setLong(1, vcf_id);
				row =  pstmt2.executeQuery();
				
				final String CHROM_prefix="#CHROM\t";
				while(row.next()) {
					final Clob clob = row.getClob(1);
					final Reader r= clob.getCharacterStream();
					/* read the first bytes to check if it's the #CHROM line
					 * if true, add a VCF header line with VCF ID and NAME
					 *  */
					final StringBuilder chrom_header_line = new StringBuilder(CHROM_prefix.length());
					int c;
					while(chrom_header_line.length()< CHROM_prefix.length() &&
						((c=r.read())!=1) )
						{
						chrom_header_line.append((char)c);
						}
					final String amorce = chrom_header_line.toString();
					if(amorce.equals(CHROM_prefix))
						{
						pwOut.println(VCF_HEADER_FILE_ID+vcf_id);
						pwOut.println(VCF_HEADER_FILE_NAME+vcfName);
						}
					pwOut.print(amorce);
					IOUtils.copyTo(r,pwOut);
					r.close();
					pwOut.println();
				}
				
				
				row.close();
				num_vcf_exported++;
				
				}
				
			pstmt.close();
			pstmt2.close();
			
			pwOut.flush();
			pwOut.close();
			
			if(num_vcf_exported==0 ) {
				LOG.warn("NO VCF WAS EXPORTED");
			}
			LOG.info("count(VCF) exported "+num_vcf_exported);
			return RETURN_OK;
		} catch (final Exception e) {
			return wrapException(e);
		} finally {
			CloserUtil.close(pwOut);
			CloserUtil.close(row);
			CloserUtil.close(pstmt);
			CloserUtil.close(pstmt2);
		}
	}

	
	
	private Collection<Throwable> doCommandDump(){
		final Pattern comma = Pattern.compile("[,]");
		final List<String> args = getInputFiles();
		final Set<Long> vcfIds = new TreeSet<>();
		try {
			for(final String token: args) {
				for(final String idStr: comma.split(token))
					{
					if(idStr.trim().isEmpty()) continue;
					LOG.info("Getting VCF "+idStr);
					long vcf_id = -1L;
					try {
						vcf_id = Long.parseLong(idStr);
					} catch(final NumberFormatException err) {
						return wrapException("Bad VCF ID :"+idStr);
					}
					vcfIds.add(vcf_id);
					}
				}
			return dump(vcfIds);
		} catch (final Exception e) {
			return wrapException(e);
		} finally {
		}
	}

	
	private Collection<Throwable> doReadConcatenatedVcf(){
		int number_of_ref_allele_truncated=0;
		PreparedStatement pstmt = null;
		PreparedStatement pstmt2 = null;
		PreparedStatement pstmt3 = null;
		ResultSet row = null;
		PrintWriter pw = null;
		final List<String> args = new ArrayList<>(IOUtils.unrollFiles(getInputFiles()));
		LOG.info(args.toString());
		LineIterator lineIter=null;
		final String titleHeaderTag = (
				super.titleHeaderStr==null || super.titleHeaderStr.trim().isEmpty()?
				null:
				"##"+titleHeaderStr+"="
				);
		try {
			int fileidx=0;
			
			pw = openFileOrStdoutAsPrintWriter();
			pw.println("#ID\tNAME");

			do
			{
				if(fileidx==0 && args.isEmpty()) {
					lineIter = IOUtils.openStreamForLineIterator(stdin());
				} else
				{
					lineIter = IOUtils.openURIForLineIterator(args.get(fileidx));
				}
				int num_vcf_in_this_stream = 0;
				while(lineIter.hasNext()) {
					String filename= "vcf"+(++ID_GENERATOR);
					if(num_vcf_in_this_stream==0 && !args.isEmpty()) {
						filename = args.get(fileidx);
					}
					
					final List<String> headerLines = new ArrayList<>();
					while(lineIter.hasNext() && lineIter.peek().startsWith("#")) {
						final String h= lineIter.next();
						if( h.startsWith(VCF_HEADER_FILE_ID) ||h.startsWith(VCF_HEADER_FILE_NAME)) {
							LOG.info("Ignoring line "+h);
							continue;
						}
						/* find filename in vcf header */
						if( titleHeaderTag!=null &&
							h.startsWith(titleHeaderTag) &&
							h.trim().length()>titleHeaderTag.length()) {
							filename = h.substring(titleHeaderTag.length()).trim();
						}
						
						headerLines.add(h);
					}
					final VCFUtils.CodecAndHeader cah = VCFUtils.parseHeader(headerLines);
					
					pstmt = this.conn.prepareStatement("INSERT INTO VCF(NAME) VALUES(?)",PreparedStatement.RETURN_GENERATED_KEYS);
					pstmt.setString(1, filename);
					if(pstmt.executeUpdate()!=1) {
						return wrapException("Cannot insert VCF ?");
					}
					final long vcf_id =getLastGeneratedId(pstmt);
					pstmt.close();
					
					pw.print(vcf_id);
					pw.print("\t");
					pw.println(filename);
					pw.flush();
					
					pstmt = this.conn.prepareStatement("SELECT ID FROM ROWCONTENT WHERE MD5SUM=?");
					pstmt2 = this.conn.prepareStatement("INSERT INTO ROWCONTENT(MD5SUM,CONTENT,CONTIG,START,STOP,ALLELE_REF,FILTERED) VALUES (?,?,?,?,?,?,?)",PreparedStatement.RETURN_GENERATED_KEYS);
					pstmt3 = this.conn.prepareStatement("INSERT INTO VCFROW(VCF_ID,ROW_ID) VALUES (?,?)");
					pstmt3.setLong(1, vcf_id);
					
					final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(cah.header);
					/* insert VCF header lines */
					for(final String line:headerLines) {
						final String md5=super.md5(line);
						long content_id = -1L;
						pstmt.setString(1, md5);
						row = pstmt.executeQuery();
						while(row.next()) {
							content_id = row.getLong(1);
						}
						row.close();				
						
						/* vcf content was not found, create it */
						if(content_id==-1L) {
							pstmt2.setString(1, md5);
							
							pstmt2.setString(2,line);
							pstmt2.setNull(3,Types.VARCHAR);
							pstmt2.setNull(4,Types.INTEGER);
							pstmt2.setNull(5,Types.INTEGER);
							pstmt2.setNull(6,Types.VARCHAR);
							pstmt2.setShort(7, (short)1);
							if(pstmt2.executeUpdate()!=1) {
								return wrapException("Cannot insert ROWCONTENT ?");
							}
							content_id =getLastGeneratedId(pstmt2);
						}
						
						/* insert new VCF row */
						pstmt3.setLong(2, content_id);
						if(pstmt3.executeUpdate()!=1) {
							return wrapException("Cannot insert VCFROW ?");
						}
					}
					
					LOG.info("Inserted "+filename+" ID="+vcf_id);
					while(lineIter.hasNext() && !lineIter.peek().startsWith("#")) {
						final String line = lineIter.next();
						final String md5 = super.md5(line);
						
						long content_id = -1L;
						pstmt.setString(1, md5);
						row = pstmt.executeQuery();
						while(row.next()) {
							content_id = row.getLong(1);
						}
						row.close();
						/* vcf variants content was not found, create it */
						if(content_id==-1L) {
							/* decode to get chrom/start/end/ref */
							final VariantContext ctx = progress.watch(cah.codec.decode(line));
							
							pstmt2.setString(1, md5);
							
							pstmt2.setString(2,line);
							pstmt2.setString(3, ctx.getContig());
							pstmt2.setInt(4, ctx.getStart());
							pstmt2.setInt(5, ctx.getEnd());
							String refBase =ctx.getReference().getBaseString();
							/* sql table for Ref_allele is a varchar(MAX_REF_BASE_LENGTH) */
							if(refBase.length()>50) {
								LOG.warn("Warning: TRUNCATING LARGE REF BASE TO FIT IN DATABASE : VARCHAR("+MAX_REF_BASE_LENGTH+") characters:"+refBase);
								refBase = refBase.substring(0,MAX_REF_BASE_LENGTH);
								++number_of_ref_allele_truncated;
							}
							pstmt2.setString(6,refBase );
							pstmt2.setShort(7, (short)(ctx.isFiltered()?1:0));
							if(pstmt2.executeUpdate()!=1) {
								return wrapException("Cannot insert ROWCONTENT ?");
							}
							content_id =getLastGeneratedId(pstmt2);
						}
						
						/* insert new VCF row */
						pstmt3.setLong(2, content_id);
						if(pstmt3.executeUpdate()!=1) {
							return wrapException("Cannot insert VCFROW ?");
						}
					}
					pstmt2.close();
					pstmt3.close();
					pstmt.close();
					progress.finish();
					num_vcf_in_this_stream++;
					} /* end of while iter has next */
				CloserUtil.close(lineIter);
				lineIter=null;
				fileidx++;
			} while(fileidx < args.size());
			
			pw.flush();
			pw.close();
			
			compress();
			LOG.warn("Number of REF alleles length(REF)> VARCHAR("+MAX_REF_BASE_LENGTH+") truncated:"+number_of_ref_allele_truncated);
			return RETURN_OK;
		} catch (final Exception e) {
			return wrapException(e);
		} finally {
			CloserUtil.close(pw);
			CloserUtil.close(row);
			CloserUtil.close(pstmt);
			CloserUtil.close(pstmt2);
			CloserUtil.close(pstmt3);
			CloserUtil.close(lineIter);
		}
	}

			
	
	private Collection<Throwable> doCommandList(){
		Statement pstmt = null;
		ResultSet row = null;
		PrintWriter pw = null;
		if(!getInputFiles().isEmpty()) {
			return wrapException("Too many arguments");
		}
		try {
			/** user can add extra columns in the VCF, we collect the names of the columns */
			final List<String> cols =new ArrayList<>();
			pstmt = this.conn.createStatement();
			row = pstmt.executeQuery("SELECT  sys.SYSCOLUMNS.COLUMNNAME FROM  sys.SYSCOLUMNS, sys.SYSTABLES WHERE  sys.SYSTABLES.TABLEID= sys.SYSCOLUMNS.REFERENCEID AND  sys.SYSTABLES.TABLENAME=\'VCF\'");
			while(row.next()) {
				cols.add("VCF."+row.getString(1));
				}
			row.close();
			pstmt.close();
			final String sql=new StringBuilder("SELECT ").
					append(String.join(",", cols)).
					append(",COUNT(VCFROW.ID)  as \"COUNT_VARIANTS\" " 
						+ "FROM VCF,VCFROW,ROWCONTENT "
						+ "WHERE VCFROW.VCF_ID=VCF.ID AND VCFROW.ROW_ID = ROWCONTENT.ID AND ROWCONTENT.CONTIG IS NOT NULL "
						+ "GROUP BY "
						).
				append(String.join(",", cols)).
				toString()
				;
			
			LOG.info(sql);
			pstmt = this.conn.createStatement();
			row = pstmt.executeQuery(sql);
			final ResultSetMetaData meta = row.getMetaData();
			pw = super.openFileOrStdoutAsPrintWriter();
			for(int i=0;i< meta.getColumnCount();i++)
				{
				pw.print((i==0?"#":"\t"));
				pw.print(meta.getColumnLabel(i+1));
				}
			pw.println();
			while(row.next()) {
				for(int i=0;i< meta.getColumnCount();i++)
				{
				pw.print((i==0?"":"\t"));
				pw.print(row.getString(i+1));
				}
			pw.println();
			}
			
			pw.flush();
			return RETURN_OK;
		} catch (Exception e) {
			return wrapException(e);
		} finally {
			CloserUtil.close(pw);
			CloserUtil.close(row);
			CloserUtil.close(pstmt);
		}
	}
	
	private Collection<Throwable> doCommandDelete(){
		final Pattern comma = Pattern.compile("[,]");
		final List<String> args = getInputFiles();
		PreparedStatement pstmt1=null;
		PreparedStatement pstmt2=null;
		try {
			pstmt1 = this.conn.prepareStatement("DELETE FROM VCFROW WHERE VCF_ID=?");
			pstmt2 = this.conn.prepareStatement("DELETE FROM VCF WHERE ID=?");
			for(final String token: args) {
				for(final String idStr: comma.split(token))
					{
					if(idStr.trim().isEmpty()) continue;
					long vcf_id=-1L;
					try {
						vcf_id = Long.parseLong(idStr);
					} catch (NumberFormatException e) {
						vcf_id=-1L;
					}
					if(vcf_id<1L) {
						LOG.warn("Bad VCF ID :"+idStr);
						continue;
					}
					pstmt1.setLong(1, vcf_id);
					pstmt2.setLong(1, vcf_id);
					LOG.info("deleting vcf id:"+vcf_id);
					pstmt1.executeUpdate();
					pstmt2.executeUpdate();
					}
				}
			compress();
			return RETURN_OK;
		} catch (final Exception e) {
			return wrapException(e);
		} finally {
			CloserUtil.close(pstmt1);
			CloserUtil.close(pstmt2);
		}
	}

	
	@Override
	public Collection<Throwable> call() throws Exception {
		try {
			openDerby();
			final String command  = String.valueOf(super.actionStr);
			
			if(command.equals("read")) {
				return doReadConcatenatedVcf();
			} else if(command.equals("list")) {
				return doCommandList();
			}else if(command.equals("dump")) {
				return doCommandDump();
			}
			else if(command.equals("dumpall")) {
				return doCommandDumpAll();
			}
			else if(command.equals("dumpuniq")) {
				return doCommandDumpUniq();
			}
			else if(command.equals("delete")) {
				return doCommandDelete();
			}
			else {
				return wrapException("unknown command : "+command);
			}			
		} catch (Exception e) {
			return wrapException(e);
		} finally {
			closeDerby();
		}
		}
	
	
	public static void main(String[] args)
		{
		new VcfDerby01().instanceMainWithExit(args);
		}
	}
