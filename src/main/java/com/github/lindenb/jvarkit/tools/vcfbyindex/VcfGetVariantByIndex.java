package com.github.lindenb.jvarkit.tools.vcfbyindex;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.OptionalLong;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFHeader;
/**
BEGIN_DOC

## Example

```bash
# get random indexes
$ gunzip -c input.vcf.gz | grep -v "#" | awk '{print NR;}' | shuf | head -n 15 > index.list
# get those 15 variants
java -jar dist/jvarkit.jar vcfgetvariantbyindex -i index.list input.vcf.gz > output.vcf

```

END_DOC
*/
@Program(name="vcfbyindex",
	description="Access a Plain or BGZF-compressed VCF file by index",
	keywords={"vcf"},
	modificationDate = "20250528",
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class VcfGetVariantByIndex extends Launcher
	{
	private static Logger LOG=Logger.of(VcfGetVariantByIndex.class);
	
	@Parameter(names="-o",description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"-i","--index-file"},description=" (file) list of 1-based indexes for query")
	private Path fileListOfIndexes=null;
	@Parameter(names={"-I","--index"},description="Comma separated of 1-based index for query")
	private String indexStr="";
	@Parameter(names={"-f","--force"},description="Force the creation of the index")
	private boolean force_index_creation=false;
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate= new WritingVariantsDelegate();

	
	
	private static final String STANDARD_EXTENSION=".ith";
	
	
	/** BASE class for index */
	private  abstract class IndexFile
		implements AutoCloseable
		{
		protected final File vcfFile;
		protected final File indexFile;
		protected VCFUtils.CodecAndHeader codecAndHeader;
		private RandomAccessFile randomAccessIndexFile;
		private long count=0L;
		public IndexFile(final File vcfFile) throws IOException
			{
			this.vcfFile=vcfFile;
			this.indexFile = new File(vcfFile.getParentFile(),vcfFile.getName()+STANDARD_EXTENSION);
			}
		public VCFHeader getHeader() {
			return codecAndHeader.header;
			}
		public AbstractVCFCodec getCodec()
			{
			return codecAndHeader.codec;
			}
		
		public abstract int read() throws IOException;
		
		protected  String readLine() throws IOException
			{
	        StringBuffer buf = new StringBuffer();
	        int c;
	        while ((c = this.read()) >= 0 && c != '\n')
	            buf.append((char) c);
	        if (c < 0) return null;
	        return buf.toString();
			}
		
		public void open() throws IOException
			{
			if(this.indexFile.exists())
				{
				openForReading();
				}
			else
				{
				if(!VcfGetVariantByIndex.this.force_index_creation) {
					throw new IOException("index by offset doesn't exists. Use option --force to enable the creation of the index "+this.indexFile);
					}
				LOG.info("building index "+this.indexFile);
				openForBuilding();
				}
			}
		
		protected abstract void openVcfFile() throws IOException;
		protected abstract long getFilePointer() throws IOException;
		
		
		private void parseHeaderAndCodec() throws IOException {
			String line=null;
			List<String> headerLines=new ArrayList<>();
			while((line=readLine())!=null)
				{
				headerLines.add(line);
				if(line.startsWith("#CHROM"))  break;
				}
			this.codecAndHeader= VCFUtils.parseHeader(headerLines);
			}
		
		private void openForBuilding() throws IOException
			{
			try
				{
				this.randomAccessIndexFile=new RandomAccessFile(this.indexFile,"rw");
				openVcfFile();
				parseHeaderAndCodec();
				this.count=0L;
				long virtualPtr=this.getFilePointer();
				while(readLine()!=null)
					{
					this.randomAccessIndexFile.writeLong(virtualPtr);
					this.count++;
					virtualPtr=this.getFilePointer();
					}
				}
			catch(Exception err)
				{
				throw new IOException(err);
				}
			}
		
		private void openForReading() throws IOException
			{
			try
				{
				this.randomAccessIndexFile = new RandomAccessFile(this.indexFile,"r");
				openVcfFile();
				parseHeaderAndCodec();
				this.count = (this.indexFile.length()/Long.BYTES);
				}
			catch(Exception err)
				{
				throw new IOException(err);
				}
			}
		
		public long getVirtualPtr(long index) throws IOException
			{
			if(index<0 || index>=this.count) throw new IndexOutOfBoundsException();
			this.randomAccessIndexFile.seek(index*Long.BYTES);
			return this.randomAccessIndexFile.readLong();
			}
		
		public abstract String getLine(long index)   throws IOException;

		
		OptionalLong parseIndex(final String line) {
			long ith;
			try {
				ith=Long.parseLong(line);
				} 
			catch (NumberFormatException e) {
				LOG.error("Bad index in "+line+" ignoring");
				return OptionalLong.empty();
				}
			ith--;//0-based index
			if(ith<0 || ith>=this.size())
				{
				LOG.error("Index out of bound in "+line+" ignoring");
				return OptionalLong.empty();
				}
			return OptionalLong.of(ith);
			}
		
		public long size()
			{
			return count;
			}
		@Override
		public void close() throws IOException {
			if(randomAccessIndexFile!=null) randomAccessIndexFile.close();
			this.randomAccessIndexFile=null;
			}
		}
	
	private  class BGZIndexFile extends IndexFile
		{
		private BlockCompressedInputStream bgzin=null;
		BGZIndexFile(File f) throws IOException
			{
			super(f);
			}
		@Override
		protected void openVcfFile() throws IOException {
			this.bgzin = new BlockCompressedInputStream(this.vcfFile);
			}
		
		@Override
		public void close() throws IOException {
			super.close();
			if(this.bgzin!=null) this.bgzin.close();
			this.bgzin=null;
			}
		@Override
		public  int read() throws IOException
			{
			return this.bgzin.read();
			}
		@Override
		public String getLine(long index)   throws IOException
			{
			long offset = getVirtualPtr(index);
			this.bgzin.seek(offset);
			return readLine();
			}
		@Override
		protected long getFilePointer()  throws IOException {
			return this.bgzin.getFilePointer();
			}
		}
	
	
	/** specialized index for plain vcf */
	private  class PlainVCFIndexFile extends IndexFile
		{
		RandomAccessFile randomeAccessVCFFile =null;
		PlainVCFIndexFile(File f) throws IOException
			{
			super(f);
			}
		
		@Override
		protected void openVcfFile() throws IOException {
			this.randomeAccessVCFFile =new RandomAccessFile(this.vcfFile, "r");
			}
		
		@Override
		public void close() throws IOException {
			super.close();
			if(this.randomeAccessVCFFile !=null) this.randomeAccessVCFFile .close();
			this.randomeAccessVCFFile =null;
			}
		@Override
		public  int read() throws IOException
			{
			return this.randomeAccessVCFFile.read();
			}
		@Override
		public String getLine(long index)   throws IOException
			{
			long offset = getVirtualPtr(index);
			this.randomeAccessVCFFile.seek(offset);
			return readLine();
			}
		@Override
		protected long getFilePointer() throws IOException {
			return this.randomeAccessVCFFile .getFilePointer();
			}

		}

	
	private IndexFile openIndex(File vcfFile) throws IOException {
		final IndexFile indexFile;
		if(vcfFile.getName().toString().endsWith(FileExtensions.COMPRESSED_VCF) || 
			vcfFile.getName().toString().endsWith("vcf.bgz") )
			{
			indexFile = new BGZIndexFile(vcfFile);
			}
		else if(vcfFile.getName().endsWith(".vcf"))
			{
			indexFile =  new PlainVCFIndexFile(vcfFile);
			}
		else
			{
			throw new IOException("Not a .vcf or .vcf.gz, .vcf.bgz file: "+vcfFile);
			}
		indexFile.open();
		return indexFile;
		}
	
	public int doWork(final List<String> args) {
		try {
			final File vcfFile= new File(oneAndOnlyOneFile(args));
			try(IndexFile indexFile = openIndex(vcfFile)) {
				try(VariantContextWriter w = this.writingVariantsDelegate.dictionary(indexFile.getHeader()).open(outputFile)) {
					final VCFHeader h2=new VCFHeader(indexFile.getHeader());
					JVarkitVersion.getInstance().addMetaData(this, h2);
					w.writeHeader(h2);
					if(fileListOfIndexes!=null) {
						try(BufferedReader r=IOUtils.openPathForBufferedReading(fileListOfIndexes)) {
							String line;
							while((line=r.readLine())!=null)
								{
								if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
								final OptionalLong ith = indexFile.parseIndex(line);
								if(ith.isEmpty()) continue;
								String varStr = indexFile.getLine(ith.getAsLong());
								w.add(indexFile.getCodec().decode(varStr));
								}
							}
						}
					if(!StringUtils.isBlank(this.indexStr)) {
						for(String stri : CharSplitter.COMMA.split(this.indexStr)) {
							if(StringUtils.isBlank(stri)) continue;
							final OptionalLong ith = indexFile.parseIndex(stri);
							if(ith.isEmpty()) continue;
							String varStr = indexFile.getLine(ith.getAsLong());
							w.add(indexFile.getCodec().decode(varStr));
							}
						}
					}
				}
			} 
		catch (Throwable e)
			{
			LOG.error(e);
			return -1;
			}
		return 0;
		}
	

	
	public static void main(String[] args) throws IOException
		{
		new VcfGetVariantByIndex().instanceMainWithExit(args);
		}
	}
