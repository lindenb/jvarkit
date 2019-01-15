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


*/
package com.github.lindenb.jvarkit.tools.vcflist;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;

import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

/**
 * Factory for offsets in vcf file
 *
 */
public class VcfOffsetsIndexFactory {
	private static final Logger LOG=Logger.build(VcfOffsetsIndexFactory.class).make();
	
	public static final String INDEX_EXTENSION =".offsets";
	static final byte MAGIC[]="vcfindex.0.1".getBytes();
	private Predicate<VariantContext> acceptVariant=null;
	private Logger logger = LOG;
	
	
	/** return true if the index file exists and we can read {@link #MAGIC} */
	public static boolean hasMagic(final File indexFile)
		{
		if(indexFile==null || !indexFile.exists() || !indexFile.canRead()) {
			return false;
			}
		try(final InputStream indexio = new FileInputStream(indexFile)) {
			final byte magic[]=new byte[VcfOffsetsIndexFactory.MAGIC.length];
			if(indexio.read(magic)!= magic.length) return false;
			return Arrays.equals(magic, VcfOffsetsIndexFactory.MAGIC);
			}
		catch(final IOException err) {
			return false;
			}
		}
	/** get default index file associated to this vcf file */
	public static File getDefaultIndexFile(final File vcf) {
		return new File(vcf.getParentFile(),vcf.getName()+INDEX_EXTENSION);
		}
	/** set a predicate for the variant that will be indexed (default is 'all' ) */
	public VcfOffsetsIndexFactory setPredicate(final Predicate<VariantContext> acceptVariant) {
		this.acceptVariant = acceptVariant;
		return this;
		}
	public VcfOffsetsIndexFactory setLogger(final Logger logger) {
		this.logger = logger;
		return this;
		}
	public File indexVcfFileIfNeeded(final File vcfFile) throws IOException
		{
		return indexVcfFileIfNeeded(vcfFile, getDefaultIndexFile(vcfFile));
		}
	/** index vcf file if needed: index file exists, contains {@link #MAGIC} and was created after vcf */
	public File indexVcfFileIfNeeded(final File vcfFile,final File indexFile) throws IOException
		{		
		if( indexFile!=null &&
			indexFile.exists() && 
			indexFile.length() >= MAGIC.length && 	
			indexFile.lastModified() >= vcfFile.lastModified() &&
			hasMagic(indexFile)) {
			return indexFile;
			}
		else
			{
			return 	indexVcfFile(vcfFile,indexFile);
			}
		}
	
	public File indexVcfFile(final File vcfFile) throws IOException
		{
		return indexVcfFile(vcfFile, getDefaultIndexFile(vcfFile));
		}
	/** index a vcf file for its variant offsets */
	public File indexVcfFile(final File vcfFile,final File indexFile) throws IOException
		{
		LOG.info("indexing "+vcfFile);
		IOUtil.assertFileIsReadable(vcfFile);
		DataOutputStream daos = null;
		BlockCompressedInputStream bgzin = null;
		AsciiLineReader ascii = null;
		VCFHeader header=null;
		
		final VCFCodec codec = new VCFCodec();
		SAMSequenceDictionaryProgress progress=null;
		
		try {
			daos = new DataOutputStream(new FileOutputStream(indexFile));
			daos.write(MAGIC);
			if(vcfFile.getName().endsWith(".vcf.gz")) {
				bgzin = new BlockCompressedInputStream(vcfFile);
				ascii = null;			
				}
			else if(vcfFile.getName().endsWith(".vcf"))
				{
				bgzin  = null;
				ascii  = new AsciiLineReader(new FileInputStream(vcfFile));
				}
			else
				{
				throw new IllegalArgumentException("not a vcf.gz or vcf file: "+vcfFile);
				}
			final List<String> headerLines=new ArrayList<>();
			for(;;)
				{
				final long offset = (ascii==null?bgzin.getPosition():ascii.getPosition());
				final String line = (ascii==null?bgzin.readLine():ascii.readLine());
				if(line==null) break;
				if(line.startsWith("#")) {
					headerLines.add(line);
					if(line.startsWith("#CHROM"))
						{
						codec.readHeader(new LineIterator() {
							int i=0;
							@Override
							public String next() {
								final String s= headerLines.get(i);
								i++;
								return s;
							}
							
							@Override
							public boolean hasNext() {
								return i < headerLines.size();
							}
							
							@Override
							public String peek() {
								return i < headerLines.size()?
										headerLines.get(i):
										null;
							}
						});
						
						header = VCFUtils.parseHeader(headerLines).header;
						progress = new SAMSequenceDictionaryProgress(header);
						progress.logger(this.logger==null?LOG:this.logger);
						progress.setLogPrefix("indexing");
						}
					continue;
					}
				if(progress==null) {
					throw new JvarkitException.FileFormatError("no vcf header in "+vcfFile);
					}
				final VariantContext ctx = codec.decode(line);
				progress.watch(ctx);
				if(this.acceptVariant!=null ) {
					if(!acceptVariant.test(ctx)) continue;
					}
				
				daos.writeLong(offset);
				}
			if(progress==null) {
				throw new JvarkitException.FileFormatError("no vcf header in "+vcfFile);
				}
			progress.finish();
			daos.flush();
			daos.close();
			return indexFile;
			}
		catch(final IOException err)
			{
			throw err;
			}
		finally
			{
			CloserUtil.close(ascii);
			CloserUtil.close(bgzin);
			CloserUtil.close(daos);
			}
		}
	}
