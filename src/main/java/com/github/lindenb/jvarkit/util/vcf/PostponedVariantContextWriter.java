package com.github.lindenb.jvarkit.util.vcf;

import java.io.File;
import java.io.PrintStream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

public class PostponedVariantContextWriter implements VariantContextWriter {
	
	public static class WritingVcfConfig
		{
		@Parameter(names={"--vcfmd5"},description="VCF, create MD5 checksum when writing a VCF/BCF to a file.")
		private boolean createMd5=false;
		@Parameter(names={"--vcfcreateindex"},description="VCF, create tribble or tabix Index when writing a VCF/BCF to a file.")
		private boolean createIndex=false;
		@Parameter(names={"--outputbcf"},description="Output bcf (for streams)")
		private boolean outputbcf=false;
	
		private SAMSequenceDictionary dict=null;

		public WritingVcfConfig(final WritingVcfConfig config) {
			this.createIndex = config.createIndex;
			this.createMd5 = config.createMd5;
			this.outputbcf = config.outputbcf;
			this.dict = config.dict;
			}
		
		public WritingVcfConfig() {
			}
		
		public WritingVcfConfig dictionary(final SAMSequenceDictionary dict)
			{
			this.dict=dict;
			return this;
			}
		
		protected VariantContextWriterBuilder createVariantContextWriterBuilder()
			{
			return new VariantContextWriterBuilder();
			}
		
		}

	private final WritingVcfConfig config;
	private final PrintStream _stdout;
	private final File outputFileOrNull;
	private VariantContextWriter delegate = null;
	public PostponedVariantContextWriter(
			final WritingVcfConfig config, 
			final PrintStream stdout,
			final File outputFileOrNull) {
		this.config = new WritingVcfConfig(config);
		this._stdout = (stdout==null?System.out:stdout);
		this.outputFileOrNull = outputFileOrNull;
		}
	
	public PostponedVariantContextWriter(
			final WritingVcfConfig config, 
			final File outputFileOrNull) {
		this(config,System.out,outputFileOrNull);
		}
	
	private VariantContextWriter getDelegate() {
		if( this.delegate == null ) {
			throw new IllegalStateException(
					"WriteHeader was never called for "+
							(this.outputFileOrNull==null?
							"<STDOUT>":
							this.outputFileOrNull.getPath()));
			}
		return this.delegate;
		}
	
	@Override
	public void writeHeader(final VCFHeader header) {
		if(this.delegate!=null) {
			throw new IllegalStateException("WriteHeader called twice");
			}
		
		SAMSequenceDictionary dict=config.dict;
		if(config.createIndex &&  dict==null) {
			dict= header.getSequenceDictionary();
			}
		
		final VariantContextWriterBuilder builder = this.config.createVariantContextWriterBuilder();
		builder.unsetOption(Options.INDEX_ON_THE_FLY);
		if(dict!=null) {
			builder.setReferenceDictionary(dict);
			}
		
		if(this.outputFileOrNull==null) {
			builder.setCreateMD5(false);
			if(this.config.outputbcf) {
				builder.setOutputBCFStream(this._stdout);
				}
			else
				{
				builder.setOutputVCFStream(this._stdout);
				}
			}
		else
			{
			builder.setCreateMD5(this.config.createMd5);
			if(this.config.createIndex)
				{
				if(dict==null|| dict.isEmpty()) {
					throw new JvarkitException.DictionaryMissing("Cannot index vcf when SamSequence dictionary missing");
					}
				builder.setOption(Options.INDEX_ON_THE_FLY);
				}
			builder.setOutputFile(this.outputFileOrNull);
			}
		this.delegate = builder.build();
		this.delegate.writeHeader(header);
		}

	@Override
	public void close() {
		CloserUtil.close(this.delegate);
	}

	@Override
	public boolean checkError() {
		return (this.delegate==null?false:this.delegate.checkError());
	}

	@Override
	public void add(final VariantContext vc) {
		getDelegate().add(vc);
	}

}
