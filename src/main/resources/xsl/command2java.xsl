<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet
	version='1.0'
	xmlns:c="http://github.com/lindenb/jvarkit/"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	>
<xsl:import href="mod.command.xsl"/>
<xsl:output method="text"/>
<xsl:param name="githash">undefined</xsl:param>
<xsl:param name="javaversion">7</xsl:param>

<xsl:template match="/">
 <xsl:apply-templates select="c:app"/>
</xsl:template>

<xsl:template match="c:app">
<xsl:apply-templates select="." mode="header"/>
package <xsl:apply-templates select="." mode="package"/>;


@javax.annotation.Generated("xslt")
public abstract class <xsl:apply-templates select="." mode="abstract-class-name"/>
	extends <xsl:choose>
		<xsl:when test="@extends"><xsl:value-of select="@extends"/></xsl:when>
		<xsl:otherwise>com.github.lindenb.jvarkit.util.command.CommandFactory</xsl:otherwise>
	</xsl:choose>
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(<xsl:apply-templates select="." mode="abstract-class-name"/>.class);
	<xsl:apply-templates select=".//c:option"/>
	
	<xsl:if test="c:output/@type='sam' or c:output/@type='bam'">
	private htsjdk.samtools.SamReader.Type outputformat= htsjdk.samtools.SamReader.Type.SAM_TYPE;
	</xsl:if>
	
	<xsl:if test="c:snippet[@id='sorting-collection'] or c:snippet[@id='tmp-dir']">
	/** list of tmp directories */
	private java.util.List&lt;java.io.File&gt; tmpDirs = null;
	
	
	/** add a temporary directory */
	public void addTmpDirectory(java.io.File dirOrFile)
		{
		if(dirOrFile==null) return;
		if(dirOrFile.isFile())
			{
			dirOrFile=dirOrFile.getParentFile();
			if(dirOrFile==null) return;
			}
		if(this.tmpDirs==null)
			{
			this.tmpDirs = new java.util.ArrayList&lt;java.io.File&gt;();
			}
		
		this.tmpDirs.add(dirOrFile);
		}
	/** returns a list of tmp directory */
	protected java.util.List&lt;java.io.File&gt; getTmpDirectories()
		{
		if(this.tmpDirs==null)
			{
			this.tmpDirs= new java.util.ArrayList&lt;java.io.File&gt;();
			}
		if(this.tmpDirs.isEmpty())
			{
			LOG.info("Adding 'java.io.tmpdir' directory to the list of tmp directories");
			this.tmpDirs.add(new java.io.File(System.getProperty("java.io.tmpdir")));
			}
		return this.tmpDirs;
		}
	</xsl:if>
		
	<xsl:if test="c:snippet[@id='sorting-collection']">
	/** When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed. */
	private int maxRecordsInRam = 500000;
	
	</xsl:if>
	
	<xsl:if test="not(@generate-output-option='false')">
		/** option outputFile */
		protected java.io.File outputFile = null;
		
		/** getter for outputFile */
		public java.io.File getOutputFile()
			{
			return this.outputFile;
			}
		
		/** setter for outputFile */
		public  void  setOutputFile( final java.io.File outputFile)
			{
			this.outputFile = outputFile;
			}
		</xsl:if>
	
	<xsl:if test="not(@generate-constructor='false')">
	/** Constructor */
	protected <xsl:apply-templates select="." mode="abstract-class-name"/>()
		{
		}
	</xsl:if>


	/** return application Name */
	@Override
	public String getName()
		{
		return "<xsl:apply-templates select="." mode="class-name"/>";
		}
	
	<xsl:choose>
		<xsl:when test="c:label">
		@Override
		public String getLabel()
			{
			return "<xsl:value-of select="c:label"/>";
			}
		</xsl:when>
		<xsl:when test="@label">
		@Override
		public String getLabel()
			{
			return "<xsl:value-of select="@label"/>";
			}
		</xsl:when>
	</xsl:choose>
	
	<xsl:choose>
		<xsl:when test="c:description">
		@Override
		public String getDescription()
			{
			return "<xsl:value-of select="c:description"/>";
			}
		</xsl:when>
		<xsl:when test="@description">
		@Override
		public String getDescription()
			{
			return "<xsl:value-of select="@description"/>";
			}
		</xsl:when>
	</xsl:choose>
	
	@Override
	protected void fillOptions(final org.apache.commons.cli.Options options)
		{
		<xsl:apply-templates select=".//c:option|.//c:options-group" mode="cli"/>
		
		<xsl:if test="not(@generate-output-option='false')">
		options.addOption(org.apache.commons.cli.Option
			.builder("o")
			.longOpt("output")
			.desc("output file. <xsl:if test="c:output/@type='sam' or c:output/@type='bam'"> extension should be .sam or .bam </xsl:if> Default: stdout")
			.argName("FILENAME")
			.hasArg(true)
			.type(org.apache.commons.cli.PatternOptionBuilder.FILE_VALUE)
			.build() );	
		</xsl:if>
		
		<xsl:if test="c:output/@type='sam' or c:output/@type='bam'">
		options.addOption(org.apache.commons.cli.Option
			.builder("formatout")
			.longOpt("formatout")
			.desc("output format : sam or bam if stdout")
			.argName("FORMAT")
			.hasArg(true)
			.type(org.apache.commons.cli.PatternOptionBuilder.STRING_VALUE)
			.build() );	
		</xsl:if>
		
		
		<xsl:if test="c:snippet[@id='sorting-collection'] or c:snippet[@id='tmp-dir']">
		options.addOption(org.apache.commons.cli.Option
			.builder("tmpdir")
			.longOpt("tmpdir")
			.desc("add tmp directory")
			.argName("DIR")
			.hasArgs()
			.type(org.apache.commons.cli.PatternOptionBuilder.FILE_VALUE)
			.build() );	
		</xsl:if>
		
		<xsl:if test="c:snippet[@id='sorting-collection']">
		options.addOption(org.apache.commons.cli.Option
			.builder("maxrecordsinram")
			.longOpt("maxrecordsinram")
			.desc("When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.")
			.argName("MAXRECORDS")
			.hasArg(true)
			.type(org.apache.commons.cli.PatternOptionBuilder.NUMBER_VALUE)
			.build() );	
		</xsl:if>
		

		
		super.fillOptions(options);
		}
	
	@Override
	protected  com.github.lindenb.jvarkit.util.command.CommandFactory.Status visit(final org.apache.commons.cli.Option opt)
		{
		<xsl:apply-templates select=".//c:option" mode="visit"/>
		<xsl:if test="not(@generate-output-option='false')">
		if(opt.getOpt().equals("o"))
			{
			java.io.File tmpf =  null;
			try { tmpf = new java.io.File(opt.getValue());}
			catch(Exception err) { LOG.error("Cannot cast "+opt.getValue()+" to output File",err); return com.github.lindenb.jvarkit.util.command.CommandFactory.Status.EXIT_FAILURE;}
			this.setOutputFile(tmpf);
			return com.github.lindenb.jvarkit.util.command.CommandFactory.Status.OK;
			}
		</xsl:if>
		<xsl:if test="c:output/@type='sam' or c:output/@type='bam'">
		if(opt.getOpt().equals("formatout"))
			{
			String formatout= opt.getValue().toLowerCase();
			if(!formatout.startsWith(".")) formatout="."+formatout;
			if( formatout.equals(".bam"))
				{
				this.outputformat = htsjdk.samtools.SamReader.Type.BAM_TYPE;
				}
			else if( formatout.equals(".sam"))
				{
				this.outputformat = htsjdk.samtools.SamReader.Type.SAM_TYPE;
				}
			else
				{
				LOG.error(formatout+" is not a valid extension.");
				return com.github.lindenb.jvarkit.util.command.CommandFactory.Status.EXIT_FAILURE;
				}
			return com.github.lindenb.jvarkit.util.command.CommandFactory.Status.OK;
			}
		</xsl:if>
		
		
		<xsl:if test="c:snippet[@id='sorting-collection'] or c:snippet[@id='tmp-dir']">
		if(opt.getOpt().equals("tmpdir"))
			{
			for(String s:opt.getValues())
				{
				try {
					java.io.File  dir = new java.io.File(s);
					this.addTmpDirectory(dir);
					}
				catch(Exception err) { LOG.error("Cannot cast "+s+" to File",err); return com.github.lindenb.jvarkit.util.command.CommandFactory.Status.EXIT_FAILURE;}
				}
			return com.github.lindenb.jvarkit.util.command.CommandFactory.Status.OK;
			}

		</xsl:if>
		
		<xsl:if test="c:snippet[@id='sorting-collection']">
		if(opt.getOpt().equals("maxrecordsinram"))
			{
			try
				{
				this.maxRecordsInRam = Integer.parseInt(opt.getValue());
				if(this.maxRecordsInRam&lt;2)
					{
					LOG.error(opt.getValue()+" is not a valid value (&lt;2).");
					}
				}
			catch(Exception err)
				{
				LOG.error(opt.getValue()+" not a valid value for "+opt.getOpt());
				return com.github.lindenb.jvarkit.util.command.CommandFactory.Status.EXIT_FAILURE;
				}
			return com.github.lindenb.jvarkit.util.command.CommandFactory.Status.OK;
			}
		</xsl:if>
		
		return super.visit(opt);
		}
		
	@Override
	public String getVersion()
		{
		return "<xsl:value-of select="$githash"/>";
		}
		
	/** Command */
	static abstract class <xsl:apply-templates select="." mode="abstract-command-name"/>
		extends
		<xsl:choose>
			<xsl:when test="@extends-command"><xsl:value-of select="@extends-command"/></xsl:when>
			<xsl:otherwise> com.github.lindenb.jvarkit.util.command.Command</xsl:otherwise>
		</xsl:choose>
		{
		protected <xsl:apply-templates select="." mode="abstract-command-name"/>()
			{
			setLog( <xsl:apply-templates select="." mode="abstract-class-name"/>.LOG);
			}
			
		@Override
		public void copyFrom(final com.github.lindenb.jvarkit.util.command.CommandFactory f) {
			<xsl:apply-templates select="." mode="abstract-class-name"/> factory = <xsl:apply-templates select="." mode="abstract-class-name"/>.class.cast(f);
			super.copyFrom(f);
			<xsl:apply-templates select=".//c:option" mode="copy"/>
			
			<xsl:if test="not(@generate-output-option='false')">
			this.setOutputFile(factory.getOutputFile());
			</xsl:if>
			<xsl:if test="c:output/@type='sam' or c:output/@type='bam'">
			this.outputformat = factory.outputformat;
			</xsl:if>
			
			<xsl:if test="c:snippet[@id='sorting-collection'] or c:snippet[@id='tmp-dir']">
			this.tmpDirs = new java.util.ArrayList&lt;java.io.File&gt;(factory.getTmpDirectories());
			</xsl:if>
		
			<xsl:if test="c:snippet[@id='sorting-collection']">
			this.maxRecordsInRam = factory.maxRecordsInRam;
			</xsl:if>
			}
			
			
		<xsl:apply-templates select=".//c:option"/>
		
		<xsl:if test="c:output/@type='sam' or c:output/@type='bam'">
		private htsjdk.samtools.SamReader.Type outputformat= htsjdk.samtools.SamReader.Type.SAM_TYPE;
		
				protected htsjdk.samtools.SAMFileWriter openSAMFileWriter(final htsjdk.samtools.SAMFileHeader header,final boolean presorted)
			{
			final htsjdk.samtools.SAMFileWriterFactory sfw= new htsjdk.samtools.SAMFileWriterFactory();
			if(getOutputFile()==null)
				{
				if(this.outputformat==null || this.outputformat.equals(htsjdk.samtools.SamReader.Type.SAM_TYPE))
					{
					LOG.info("Saving as SAM");
					return sfw.makeSAMWriter(header, presorted, stdout());
					}
				else if( this.outputformat.equals(htsjdk.samtools.SamReader.Type.BAM_TYPE))
					{
					LOG.info("Saving as BAM");
					return sfw.makeBAMWriter(header, presorted, stdout());
					}
				else
					{
					throw new IllegalStateException("Bad output format");
					}
				}
			else
				{
				LOG.info("Saving as "+ getOutputFile());
				return sfw.makeSAMOrBAMWriter(header, presorted, getOutputFile());
				}
			}

		
		</xsl:if>
		
		
		
		
		<xsl:if test="not(@generate-output-option='false')">
		/** option outputFile */
		protected java.io.File outputFile = null;
		
		/** getter for outputFile */
		public java.io.File getOutputFile()
			{
			return this.outputFile;
			}
		
		/** setter for outputFile */
		public  void  setOutputFile( final java.io.File outputFile)
			{
			this.outputFile = outputFile;
			}
		
		protected java.io.PrintStream openFileOrStdoutAsPrintStream() throws java.io.IOException
			{
			if(getOutputFile()!=null)
				{
				return new java.io.PrintStream(getOutputFile());
				}
			else
				{
				return stdout();
				}
			}
		protected java.io.OutputStream openFileOrStdoutAsStream() throws java.io.IOException
			{
			if(getOutputFile()!=null)
				{
				return com.github.lindenb.jvarkit.io.IOUtils.openFileForWriting(getOutputFile());
				}
			else
				{
				return stdout();
				}
			}
		</xsl:if>
		
		@Override
		public void cleanup()
			{
			super.cleanup();
			}
		
		<xsl:choose>
		<xsl:when test="c:output/@type='vcf'">
		/* creates a VariantContextWriter according to FileOUt */
		protected htsjdk.variant.variantcontext.writer.VariantContextWriter  openVariantContextWriter()
			throws java.io.IOException
			{
			if(getOutputFile()!=null)
				{
				return com.github.lindenb.jvarkit.util.vcf.VCFUtils.createVariantContextWriter(this.getOutputFile());
				}
			else
				{
				return com.github.lindenb.jvarkit.util.vcf.VCFUtils.createVariantContextWriterToOutputStream(stdout());
				}
			}
		</xsl:when>
		</xsl:choose>
		

		
		<xsl:choose>
		<xsl:when test="not(c:input/@type)">
			<xsl:message terminate="no">warning: input type undefined</xsl:message>
		</xsl:when>
		<xsl:when test="c:input/@type='stdin-or-one' or c:input/@type='sam' or c:input/@type='vcf'">
		

		<xsl:if test="c:input/@type='vcf'">
		
		/* creates a VCF iterator from inputName. */
		protected com.github.lindenb.jvarkit.util.vcf.VcfIterator openVcfIterator(final String inputName)
			throws java.io.IOException
			{
			if( inputName == null )
					{
					return com.github.lindenb.jvarkit.util.vcf.VCFUtils.createVcfIteratorFromStream(stdin());
					}
				else
					{
					return com.github.lindenb.jvarkit.util.vcf.VCFUtils.createVcfIterator(inputName);
					}
			}

		</xsl:if>
		
		<xsl:if test="c:input/@type='sam'">
		
		protected htsjdk.samtools.SamReaderFactory createSamReaderFactory()
			{
			return  htsjdk.samtools.SamReaderFactory.makeDefault().validationStringency(htsjdk.samtools.ValidationStringency.LENIENT);
			}
		
		protected htsjdk.samtools.SamReader openSamReader(final String inputName)
			{
			final htsjdk.samtools.SamReaderFactory srf= this.createSamReaderFactory();
			if(inputName==null)
				{
				LOG.info("opening stdin");
				return srf.open(htsjdk.samtools.SamInputResource.of(stdin()));
				}
			else
				{
				LOG.info("opening "+inputName);
				return srf.open(htsjdk.samtools.SamInputResource.of(inputName));
				}
			}

		</xsl:if>
		
		/** program should process this file or stdin() if inputName is null */ 
		protected abstract java.util.Collection&lt;Throwable&gt; call(final String inputName) throws Exception;
		
		@Override
		public  java.util.Collection&lt;Throwable&gt; call() throws Exception
			{
			final java.util.List&lt;String&gt; args= getInputFiles();
			if(args.isEmpty())
				{
				return call(null);
				}
			else if(args.size()==1)
				{
				final String filename = args.get(0);
				return call(filename);
				}
			else
				{
				return wrapException(getMessageBundle("illegal.number.of.arguments"));
				}
			}
		</xsl:when>
		<xsl:when test="c:input/@type='strings'">
		/* input type is 'strings' */
		</xsl:when>
		
		<xsl:otherwise>
		<xsl:message terminate="yes">undefined input/@type </xsl:message>
		</xsl:otherwise>
		</xsl:choose>
		
		
		<xsl:if test="c:snippet[@id='sorting-collection'] or c:snippet[@id='tmp-dir']">
		/** list of tmp directories */
		protected java.util.List&lt;java.io.File&gt; tmpDirs = new java.util.ArrayList&lt;java.io.File&gt;();
		</xsl:if>
		
		<xsl:if test="c:snippet[@id='sorting-collection']">
		/** When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed. */
		protected int maxRecordsInRam = 500000;
		</xsl:if>
		
		}
	<xsl:if test="number($javaversion) &gt;= 8">
	<xsl:apply-templates select="." mode="jfx"/>
	</xsl:if>
	
	
	}
</xsl:template>




</xsl:stylesheet>


