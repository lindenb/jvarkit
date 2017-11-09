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
package com.github.lindenb.jvarkit.tools.vcflist;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
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
	
	public static File getIndexFile(final File vcf) {
		return new File(vcf.getParentFile(),vcf.getName()+INDEX_EXTENSION);
		}
	public void setPredicate(final Predicate<VariantContext> acceptVariant) {
		this.acceptVariant = acceptVariant;
		}
	public void setLogger(final Logger logger) {
		this.logger = logger;
		}
	
	public File indexVcfFileIfNeeded(final File vcfFile) throws IOException
		{
		final File indexFile = VcfOffsetsIndexFactory.getIndexFile(vcfFile);
		
		if( indexFile.exists() && 
			indexFile.length() >= MAGIC.length && 	
			indexFile.lastModified() >= vcfFile.lastModified()) {
			return indexFile;
			}
		else
			{
			return 	indexVcfFile(vcfFile);
			}
		}
	public File indexVcfFile(final File vcfFile) throws IOException
		{
		IOUtil.assertFileIsReadable(vcfFile);
		final File indexFile = VcfOffsetsIndexFactory.getIndexFile(vcfFile);
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
				throw new IOException("not a vcf.gz or vcf file: "+vcfFile);
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
			indexFile.delete();
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
