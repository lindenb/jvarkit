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
package com.github.lindenb.jvarkit.tools.vcfviewgui;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.nio.charset.Charset;

import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public abstract class VcfFile implements NgsFile<VCFHeader,VariantContext>{
    private static final Logger LOG= Logger.build(VcfFile.class).make();
    private PedFile pedigree;
    
	protected VcfFile(final PedFile pedigree)  throws IOException
		{
		this.pedigree=(pedigree==null?PedFile.getEmptyInstance():pedigree);
		}
	
	@Override
	public abstract VcfFile reOpen() throws IOException;
	
	@Override
	public SAMSequenceDictionary getSequenceDictionary() {
		return getHeader().getSequenceDictionary();
	}

	/** get the associated pedigree, can be empty but never null */
	public PedFile getPedigree()
		{
		return this.pedigree;
		}
	
	public static VcfFile newInstance(final String s) throws IOException {
		return IOUtil.isUrl(s)?
				newInstance(new URL(s)):
				newInstance(new File(s))
				;
	}
	
	public static VcfFile newInstance(final File f) throws IOException {
		IOUtil.assertFileIsReadable(f);	
		
		PedFile ped=PedFile.getEmptyInstance();
		String pedName=null;
		if(f.getName().endsWith(".vcf.gz"))
			{
			pedName=f.getName().substring(0,f.getName().length()-7)+ PedFile.EXTENSION;
			}
		else if(f.getName().endsWith(".vcf"))
			{
			pedName=f.getName().substring(0,f.getName().length()-4)+ PedFile.EXTENSION;
			}
		final File pedFile=new File(f.getParentFile(),pedName);
		if(pedFile.exists() && pedFile.isFile() && pedFile.canRead())
			{
			try {
				final PedFile p= PedFile.load(pedFile);
				ped=p;
				}
			catch (final IOException err) {
				ped=PedFile.getEmptyInstance();
				}
			finally
				{
				}
			}
		
		return new LocalVcf(f,ped);
		}
	
	public static VcfFile newInstance(final URL url) throws IOException {		
		final File tbiFile=File.createTempFile("tmp.", ".vcf.gz.tbi");
    	tbiFile.deleteOnExit();
    	final String tabixurl=url.toExternalForm()+".tbi";
    	InputStream in=null;
    	FileOutputStream out=null;
		try {
			LOG.info("trying "+tabixurl);
			in = new URL(tabixurl).openStream();
			out = new FileOutputStream(tbiFile);
			IOUtil.copyStream(in, out);
			out.flush();
			out.close();
			in.close();
		} catch (final IOException err) {
			tbiFile.delete();
			throw new IOException("Cannot fetch "+tabixurl,err);
			}
		finally
			{
			CloserUtil.close(out);
			CloserUtil.close(in);
			}
		PedFile ped=PedFile.getEmptyInstance();
		if(url.toExternalForm().endsWith(".vcf.gz"))
			{
			String pedurl= url.toExternalForm().substring(0,url.toExternalForm().length()-7)+ PedFile.EXTENSION;
			LOG.info("trying "+pedurl);
			try {
				in = new URL( pedurl).openStream();
				final PedFile pedfile= PedFile.load(new BufferedReader(new InputStreamReader(in, Charset.forName("UTF-8"))));
				ped=pedfile;
				}
			catch (final IOException err) {
				ped=PedFile.getEmptyInstance();
				}
			finally
				{
				CloserUtil.close(in);
				}
			}
		return new RemoteVcfFile(url,tbiFile,ped);
		}
	
	
	private static class LocalVcf extends VcfFile {
	private final VCFFileReader vcfFileReader;
	private final File file;
	LocalVcf(final File f,final PedFile ped) throws IOException
		{
		super(ped);
		this.file=f;
		this.vcfFileReader = new VCFFileReader(f, true);
		}
	
	@Override
	public VCFHeader getHeader() {
		return this.vcfFileReader.getFileHeader();
	}


	@Override
	public CloseableIterator<VariantContext> iterator() throws IOException {
		return this.vcfFileReader.iterator();
	}

	@Override
	public CloseableIterator<VariantContext> iterator(final String contig, int start, int end) throws IOException {
		return this.vcfFileReader.query(contig,start,end);
	}

	public void close() {
		this.vcfFileReader.close();
		}
	
	@Override
		public String getSource() {
			return this.file.getAbsolutePath();
		}
	
	@Override
	public VcfFile reOpen() throws IOException {
			return new LocalVcf(this.file,getPedigree());
		}
	
	}
	
	
	private static class RemoteVcfFile extends VcfFile {
	private boolean delete_index_on_close=true;
	private final FeatureReader<VariantContext> reader;
	private final File indexFile;
	private final URL url;
	RemoteVcfFile(final URL url,final File indexFile,final  PedFile pedfile) throws IOException
		{
		super(pedfile);
		this.url = url;
		this.reader=AbstractFeatureReader.
				getFeatureReader(
						url.toExternalForm(),
						indexFile.getAbsolutePath(),
						new VCFCodec(),
						true
						)
				;
		this.indexFile=indexFile;
		}
		
	@Override
	public String getSource() {
		return url.toExternalForm();
	}
	
	@Override
		public VCFHeader getHeader() {
			return VCFHeader.class.cast(this.reader.getHeader());
			}

	@Override
	public CloseableIterator<VariantContext> iterator(String contig, final int start, int end) throws IOException {
		return this.reader.query(contig, start, end);
		}
	
	@Override
	public CloseableIterator<VariantContext> iterator() throws IOException {
		return this.reader.iterator();
	}

	@Override
	public void close() {
		if(this.delete_index_on_close)
			{
			this.indexFile.delete();	
			}
		CloserUtil.close(this.reader);
		}

	@Override
		public VcfFile reOpen() throws IOException {
			final RemoteVcfFile bf= new RemoteVcfFile(this.url, this.indexFile,getPedigree());
			bf.delete_index_on_close=false;
			return bf;
			}
	}
}
