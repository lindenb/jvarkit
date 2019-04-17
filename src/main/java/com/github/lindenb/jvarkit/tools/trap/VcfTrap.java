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
package com.github.lindenb.jvarkit.tools.trap;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Algorithms;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import htsjdk.variant.vcf.VCFIterator;

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
keywords={"vcf","trap","annotation"},
modificationDate="20190308"
)
public class VcfTrap extends Launcher {
	private static final Logger LOG = Logger.build(VcfTrap.class).make();
	@Parameter(names={"-o","--out"},required=false,description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	
	@Parameter(names={"-m","--manifest"},description="Manifest file. A tab delimited file with two columns : chromosome(tab)path-to-file-indexed-with-trapindex.",required=true)
	private Path manifestFile=null;
	@Parameter(names={"--ignore-filtered"},description="Ignore FILTERed variants (faster)")
	private boolean ignore_filtered=false;
	@Parameter(names={"-A","--attribute"},description="VCF INFO attribute Format:(ALT|GENE|SCORE)")
	private String ATT="TRAP";
	
	private static class IndexFile extends AbstractList<TrapRecord>
		implements Closeable
		{
		final String contig;
		final Path file;
		private RandomAccessFile io;
		final int _size;
		IndexFile(final String contig,final Path file) throws IOException {
			this.contig = contig;
			this.file=file;
			long length = Files.size(file);
			
			if(length< TrapIndexer.MAGIC.length)
				{
				throw new RuntimeIOException("size lower than TrapIndexer.MAGIC "+file);
				}
			length -=  TrapIndexer.MAGIC.length;
			
			if(length % TrapIndexer.RECORD_SIZOF!=0) throw new  IOException("not a multiple of "+TrapIndexer.RECORD_SIZOF+":"+length);
			this._size = (int)(length/TrapIndexer.RECORD_SIZOF);
			this.io = new RandomAccessFile(this.file.toFile(), "r");
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
	
	
	public VcfTrap() {
		}
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator iter,
			final VariantContextWriter out
			) {
		final Map<String, Path> chromToFile= new HashMap<>();

		IOUtil.assertFileIsReadable(this.manifestFile);
		BufferedReader r=null;
		try {
			r=IOUtils.openPathForBufferedReading(this.manifestFile);
			r.lines().filter(L->!L.trim().isEmpty()).
				filter(L->!L.startsWith("#")).
				forEach(L->{
				final int tab  = L.indexOf('\t');
				if(tab==-1) throw new RuntimeIOException("tab missing in "+L);
				final String contig = L.substring(0,tab);
				if(chromToFile.containsKey(contig)) {
					throw new RuntimeIOException("duplicate contig "+L);
					}
				final Path indexFile= Paths.get(L.substring(tab+1));
				IOUtil.assertFileIsReadable(indexFile);
				chromToFile.put(contig,indexFile);
				});
			}
		catch(final IOException err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			}
		
		final Function<String,IndexFile> getIndexFile = s->
			{
			Path file = chromToFile.get(s);
			if(file==null && s.startsWith("chr")) file =  chromToFile.get(s.substring(3));
			if(file==null && !s.startsWith("chr")) file = chromToFile.get("chr"+s);
			if(file==null) return null;
			try {
				return new IndexFile(s,file);
			} catch (final IOException err) {
				throw new RuntimeIOException(err);
				}
			};
		
		IndexFile current=null;
		final String ATT_MIN = this.ATT+"_MIN";
		final String ATT_MAX = this.ATT+"_MAX";
		final Set<String> contigs_not_found=new HashSet<>();
		final Comparator<TrapRecord> comparator = (A,B) ->{
			if(! A.getContig().equals(B.getContig())) throw new IllegalStateException("not the same contigs ???");
			return Integer.compare(A.getStart(), B.getStart());
			};
		
		final VCFHeader header=new VCFHeader(iter.getHeader());
	
		final VCFHeader header2=new VCFHeader(header);
		header2.addMetaDataLine(
				new VCFInfoHeaderLine(this.ATT,
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"TRAP Score:((ALT|GENE|SCORE) in Trap Database  http://trap-score.org/"
					));
		header2.addMetaDataLine(
				new VCFInfoHeaderLine(ATT_MIN,
						1,
						VCFHeaderLineType.Float,
						"Min Score in Trap Database  http://trap-score.org/"
					));
		header2.addMetaDataLine(
				new VCFInfoHeaderLine(ATT_MAX,
						1,
						VCFHeaderLineType.Float,
						"Max Score in Trap Database  http://trap-score.org/"
					));
		JVarkitVersion.getInstance().addMetaData(this, header2);
		final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().
							dictionary(iter.getHeader()).
							logger(LOG).
							build();
		out.writeHeader(header2);
		while(iter.hasNext())
			{
			final VariantContext var = progress.apply(iter.next());
			
			
			if(this.ignore_filtered && var.isFiltered())
			{
			out.add(var);
			continue;
			}
		if(current==null || !current.contig.equals(var.getContig()))
			{
			if(current!=null) 
				{
				CloserUtil.close(current);
				current=null;
				}
			if(contigs_not_found.contains(var.getContig()))
				{
				out.add(var);
				continue;
				}
			current = getIndexFile.apply(var.getContig());
			}
		if(current==null)
			{
			LOG.warn("Not indexed in trap "+var.getContig());
			contigs_not_found.add(var.getContig());
			out.add(var);
			continue;
			}
		final Set<String> annotations=new HashSet<>();
		final Float min_score[]=new Float[] {null};
		final Float max_score[]=new Float[] {null};
		
		Algorithms.equal_range_stream(
				current,
				0,
				current.size(),
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
				comparator,
				Function.identity()
				).
				filter(R->var.getReference().equals(Allele.create((byte)R.getRef(),true))).
				filter(R->var.getAlternateAlleles().stream().anyMatch(A->A.equals(Allele.create((byte)R.getAlt(),false)))).
				forEach(R->{
					
					annotations.add(String.join("|",
							String.valueOf(R.getAlt()),
							R.getGene(),
							String.format("%."+TrapIndexer.SCORE_STRLEN+"f", R.getScore())
							));
					if(min_score[0]==null || min_score[0].compareTo(R.getScore())>0)
						{
						min_score[0]=R.getScore();
						}
					if(max_score[0]==null || max_score[0].compareTo(R.getScore())<0)
						{
						max_score[0]=R.getScore();
						}
					});
			if(annotations.isEmpty())
				{
				out.add(var);
				continue;
				}
			final VariantContextBuilder vcb = new VariantContextBuilder(var);
			vcb.attribute(this.ATT, new ArrayList<>(annotations));
			vcb.attribute(ATT_MIN,min_score[0]);
			vcb.attribute(ATT_MAX,max_score[0]);
			out.add(var);
			}
		out.close();
		progress.close();
		CloserUtil.close(current);
		return 0;
		}
	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			return doVcfToVcf(args, this.outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		
		}
		
	public static void main(final String[] args) {
		new VcfTrap().instanceMainWithExit(args);
	}
}
