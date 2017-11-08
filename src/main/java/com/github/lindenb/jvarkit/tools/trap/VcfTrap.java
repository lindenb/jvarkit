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
package com.github.lindenb.jvarkit.tools.trap;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;
import javax.xml.bind.annotation.XmlType;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Algorithms;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.PostponedVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
/**
BEGIN_DOC

## Example:

```
java -jar dist/trapindexer.jar  -o chr22.dat chr22.TraPv2.txt.gz
echo -e "22\tchr22.dat" > out.manifest
java -jar dist/vcftrap.jar -m out.manifest input.vcf
```


## See also

* TrapIndexer

END_DOC
 */
@Program(name="vcftrap",
description="annotate vcf with trap database http://trap-score.org/",
keywords={"vcf","trap","annotation"}
)
public class VcfTrap extends Launcher {
	private static final Logger LOG = Logger.build(VcfTrap.class).make();
	@Parameter(names={"-o","--out"},required=false,description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@ParametersDelegate
	private PostponedVariantContextWriter.WritingVcfConfig writingVcfArgs = new PostponedVariantContextWriter.WritingVcfConfig();
	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();
	
	
	private static class IndexFile extends AbstractList<TrapRecord>
		implements Closeable
		{
		final String contig;
		final File file;
		private RandomAccessFile io;
		final int _size;
		IndexFile(final String contig,final File file) throws IOException {
			this.contig = contig;
			this.file=file;
			long length = file.length();
			
			if(length< TrapIndexer.MAGIC.length)
				{
				throw new RuntimeIOException("size lower than TrapIndexer.MAGIC "+file);
				}
			length -=  TrapIndexer.MAGIC.length;
			
			if(length % TrapIndexer.RECORD_SIZOF!=0) throw new  IOException("not a multiple of "+TrapIndexer.RECORD_SIZOF+":"+length);
			this._size = (int)(length/TrapIndexer.RECORD_SIZOF);
			this.io = new RandomAccessFile(this.file, "r");
			final byte magic[]=new byte[TrapIndexer.MAGIC.length];
			this.io.readFully(magic);
			if(!Arrays.equals(magic,  TrapIndexer.MAGIC))
				{
				this.io.close();
				throw new IOException("not a TrapIndexer file:"+file);
				}
			}
		@Override
		public TrapRecord get(final int index) {
			try {
				final byte array[]=new byte[TrapIndexer.RECORD_SIZOF];
				this.io.seek((long)TrapIndexer.MAGIC.length + (long)TrapIndexer.RECORD_SIZOF*(long)index);
				this.io.readFully(array);
				return TrapIndexer.decode(this.contig, array);
				}
			catch(final IOException err)
				{
				throw new RuntimeIOException(err);
				}
			}
		
		@Override
		public int size() {
			return this._size;
			}
		
		@Override
		public void close() throws IOException {
			LOG.debug("closing "+contig);
			CloserUtil.close(this.io);
			this.io = null;
			}
		}
	
	@XmlType(name="vcftrap")
	@XmlRootElement(name="vcftrap")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class CtxWriterFactory 
		implements VariantContextWriterFactory
			{

			@XmlElement(name="manifest")
			@Parameter(names={"-m","--manifest"},description="Manifest file. A tab delimited file with two columns : chromosome(tab)path-to-file-indexed-with-trapindex.",required=true)
			private File manifestFile=null;
			@Parameter(names={"--ignore-filtered"},description="Ignore FILTERed variants (faster)")
			private boolean ignore_filtered=false;
			@XmlElement(name="attribute")
			@Parameter(names={"-A","--attribute"},description="VCF INFO attribute Format:(ALT|GENE|SCORE)")
			private String ATT="TRAP";

			
			@XmlTransient
			final Map<String, File> chromToFile= new HashMap<>();
			
			private class CtxWriter extends DelegateVariantContextWriter
				{					
				private IndexFile current=null;
				private final boolean ignore_filtered = CtxWriterFactory.this.ignore_filtered;
				private final String ATT = CtxWriterFactory.this.ATT;
				private final String ATT_BEST = CtxWriterFactory.this.ATT+"_BEST";
				private final Set<String> contigs_not_found=new HashSet<>();
				private final Comparator<TrapRecord> comparator = (A,B) ->{
					if(! A.getContig().equals(B.getContig())) throw new IllegalStateException("not the same contigs ???");
					return Integer.compare(A.getStart(), B.getStart());
					};
				
				CtxWriter(final VariantContextWriter delegate) {
					super(delegate);
					}
				@Override
				public void writeHeader(final VCFHeader header) {
					final VCFHeader header2=new VCFHeader(header);
					header2.addMetaDataLine(
							new VCFInfoHeaderLine(this.ATT,
									VCFHeaderLineCount.UNBOUNDED,
									VCFHeaderLineType.String,
									"TRAP Score:((ALT|GENE|SCORE) in Trap Database  http://trap-score.org/"
								));
					header2.addMetaDataLine(
							new VCFInfoHeaderLine(this.ATT_BEST,
									1,
									VCFHeaderLineType.Float,
									"Best Score in Trap Database  http://trap-score.org/"
								));
					super.writeHeader(header2);
					}
				@Override
				public void add(final VariantContext var) {
					if(this.ignore_filtered && var.isFiltered())
						{
						super.add(var);
						return;
						}
					if(this.current==null || !this.current.contig.equals(var.getContig()))
						{
						if(this.current!=null) 
							{
							CloserUtil.close(this.current);
							this.current=null;
							}
						if(this.contigs_not_found.contains(var.getContig()))
							{
							super.add(var);
							return;
							}
						this.current = CtxWriterFactory.this.getIndexFile(var.getContig());
						}
					if(this.current==null)
						{
						LOG.warn("Not indexed in trap "+var.getContig());
						this.contigs_not_found.add(var.getContig());
						super.add(var);
						return;
						}
					final Set<String> annotations=new HashSet<>();
					final Float best[]=new Float[] {null};
					
					Algorithms.equal_range_stream(
							this.current,
							0,
							this.current.size(),
							new TrapRecord() {
								@Override
								public int getStart() { return var.getStart(); }
								@Override
								public int getEnd() { return var.getEnd(); }
								@Override
								public String getContig() { return var.getContig();}
								@Override
								public String getChr() { return getContig(); }
								@Override
								public float getScore() { return 0f; }
								@Override
								public char getRef() {return '\0';}
								@Override
								public String getGene() {return "";}
								@Override
								public char getAlt() { return '\0'; }
								},
							this.comparator
							).
							filter(R->var.getReference().equals(Allele.create((byte)R.getRef(),true))).
							filter(R->var.getAlternateAlleles().stream().anyMatch(A->A.equals(Allele.create((byte)R.getAlt(),false)))).
							forEach(R->{
								
								annotations.add(String.join("|",
										String.valueOf(R.getAlt()),
										R.getGene(),
										String.format("%."+TrapIndexer.SCORE_STRLEN+"f", R.getScore())
										));
								if(best[0]==null || best[0].compareTo(R.getScore())<0)
									{
									best[0]=R.getScore();
									}
								});
					if(annotations.isEmpty())
						{
						super.add(var);
						return;
						}
					final VariantContextBuilder vcb = new VariantContextBuilder(var);
					vcb.attribute(this.ATT, new ArrayList<>(annotations));
					vcb.attribute(this.ATT_BEST,best[0]);
					super.add(vcb.make());
					}
				@Override
				public void close() {
					try { this.current.close();}
					catch(IOException err) {LOG.error(err);} //CloserUtil doenst' work ??
					this.current=null;
					super.close();
					}
				}
			
			private IndexFile getIndexFile(final String s)
				{
				File file = this.chromToFile.get(s);
				if(file==null && s.startsWith("chr")) file =  this.chromToFile.get(s.substring(3));
				if(file==null && !s.startsWith("chr")) file =  this.chromToFile.get("chr"+s);
				if(file==null) return null;
				try {
					return new IndexFile(s,file);
				} catch (final IOException err) {
					throw new RuntimeIOException(err);
					}
				}
			
			@Override
			public int initialize() {
				IOUtil.assertFileIsReadable(this.manifestFile);
				BufferedReader r=null;
				try {
					r=IOUtils.openFileForBufferedReading(this.manifestFile);
					r.lines().filter(L->!L.trim().isEmpty()).
						filter(L->!L.startsWith("#")).
						forEach(L->{
						final int tab  = L.indexOf('\t');
						if(tab==-1) throw new RuntimeIOException("tab missing in "+L);
						final String contig = L.substring(0,tab);
						if(this.chromToFile.containsKey(contig)) {
							throw new RuntimeIOException("duplicate contig "+L);
							}
						final File indexFile= new File(L.substring(tab+1));
						IOUtil.assertFileIsReadable(indexFile);
						this.chromToFile.put(contig,indexFile);
						});
					}
				catch(IOException err) {
					LOG.error(err);
					return -1;
					}
				finally
					{
					CloserUtil.close(r);
					}
				return 0;
				}
			
			@Override
			public VariantContextWriter open(final VariantContextWriter delegate) {
				return new CtxWriter(delegate);
				}
			@Override
			public void close() throws IOException {
				this.chromToFile.clear();
				}
			}

	
	public VcfTrap() {
		}
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VcfIterator iter,
			final VariantContextWriter delegate
			) {	
		final VariantContextWriter out = this.component.open(delegate);
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(iter.getHeader()).logger(LOG);
		out.writeHeader(iter.getHeader());
		while(iter.hasNext())
			{
			out.add(progress.watch(iter.next()));
			}
		out.close();
		progress.finish();
		return 0;
		}
	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			if(this.component.initialize()!=0) return -1;
			return doVcfToVcf(args, this.outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.component);
			}
		}
		
	public static void main(final String[] args) {
		new VcfTrap().instanceMainWithExit(args);
	}
}
