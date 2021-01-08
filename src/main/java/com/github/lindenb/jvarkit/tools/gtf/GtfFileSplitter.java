/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.gtf;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.SortingCollection;

/**

BEGIN_DOC

### Example

split per gene

```
$ java -jar dist/gtfsplitter.jar -m gene -o TMP src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz
$ find TMP/ -name "*.gtf"
TMP/51/4db9b1fbf7674369310d8698c987b6/22_ENSG00000100403.gtf
TMP/b1/315320e53be12e62a70acb6f26a2a6/3_ENSG00000183873.gtf
TMP/b1/035d795b2e064aa90cde9c64d17014/1_ENSG00000134250.gtf

$ cat TMP/b1/315320e53be12e62a70acb6f26a2a6/3_ENSG00000183873.gtf | cut -f 3 | sort | uniq -c 
    266 CDS
    282 exon
     17 five_prime_utr
      1 gene
     11 start_codon
     10 stop_codon
      9 three_prime_utr
     14 transcript

# check same number in the GTF file for this gene

$ gunzip -c src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz | grep ENSG00000183873  | cut -f 3 | sort | uniq -c
    266 CDS
    282 exon
     17 five_prime_utr
      1 gene
     11 start_codon
     10 stop_codon
      9 three_prime_utr
     14 transcript
 
```

split per transcript:

```
$ java -jar dist/gtfsplitter.jar -m transcript -o TMP src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz
$ find TMP/ -name "*.gtf"

TMP/7f/fd225f553e6fae72d87a26d21db3e6/3_ENST00000464652.gtf
TMP/ee/4b029d630d9d04ea90c5b6d58d6476/3_ENST00000455624.gtf
(...)
TMP/d9/e65eed61df7e80fdc9a12301b3b37d/22_ENST00000352645.gtf
TMP/c9/214537e6b2363128029fdd7898d1c7/1_ENST00000479412.gtf
TMP/23/750b00f3d4b5760ebaeaeda9e39e7f/3_ENST00000327956.gtf


$ cut -f 3 TMP/7f/fd225f553e6fae72d87a26d21db3e6/3_ENST00000464652.gtf
gene
exon
transcript
exon

# check the original file (the gene is ignored with this simple grep )
$ gunzip -c src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz | grep ENST00000464652 | cut -f 3 
exon
transcript
exon

```

split per chromosome/contig

```
$ java -jar dist/gtfsplitter.jar -m contig -o TMP src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz
$ find TMP/ -name "*.gtf"
TMP/c4/ca4238a0b923820dcc509a6f75849b/1.gtf
TMP/b6/d767d2f8ed5d21a44b0e5886680cb9/22.gtf
TMP/ec/cbc87e4b5ce2fe28308fd9f2a7baf3/3.gtf

$ cut -f 1 TMP/b6/d767d2f8ed5d21a44b0e5886680cb9/22.gtf | uniq
22

```

per group

```
$ gunzip -c src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz | grep -v "#" | cut -f 3 | grep gene | wc
3

$ java -jar dist/gtfsplitter.jar -m group -p 2 -o TMP src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz
$ find TMP/ -name "*.gtf" -exec grep -w -c -H gene '{}' ';'
TMP/6f/ab500a2f6fee4cce61739dd97b11cb/Group000000.gtf:1
TMP/47/ee0407c2da530f87841e134fa6717a/Group000001.gtf:2

# with a larger gtf file

$ java -jar dist/gtfsplitter.jar -m group -p 10 -o TMP ~/jeter.gtf.gz  
INFO	2019-08-07 11:14:00	SortingCollection	Creating merging iterator from 58 files

$ find TMP/ -name "*.gtf"    -exec grep -w -c gene -H '{}' ';'
TMP/ad/7732b2eb185c413c329c1f124bb809/Group000006.gtf:6229
TMP/6f/ab500a2f6fee4cce61739dd97b11cb/Group000000.gtf:6229
TMP/38/b153383e90aaa85a260d3560d36725/Group000005.gtf:6229
TMP/4b/6c5b9a6d1830f91089bdebe8321e90/Group000002.gtf:6230
TMP/2d/e7f75cf324bb4968ef397601e67fa7/Group000003.gtf:6229
TMP/93/2ea97d325f88f3c572d8af057b1c06/Group000008.gtf:6229
TMP/96/27eea1b15190b6a726b98af9f27b3a/Group000004.gtf:6229
TMP/96/ab60170db743f9704390147fadd705/Group000007.gtf:6229
TMP/47/ee0407c2da530f87841e134fa6717a/Group000001.gtf:6230
TMP/45/a4780078f03085345ac98b5b2a766e/Group000009.gtf:6229

```

per stack

```
$ java -jar dist/gtfsplitter.jar -m stack -p 3 -o TMP src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz
$ find TMP/ -name "*.gtf" 
TMP/6f/ab500a2f6fee4cce61739dd97b11cb/Group000000.gtf
lindenb@mcclintock:~/src/jvarkit-git$ grep -w gene -c TMP/6f/ab500a2f6fee4cce61739dd97b11cb/Group000000.gtf
3

$ java -jar dist/gtfsplitter.jar -m stack -p 1 -o TMP src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz
$ find TMP/ -name "*.gtf" 
TMP/6f/ab500a2f6fee4cce61739dd97b11cb/Group000000.gtf
TMP/4b/6c5b9a6d1830f91089bdebe8321e90/Group000002.gtf
TMP/47/ee0407c2da530f87841e134fa6717a/Group000001.gtf

# with a larger gtf file

$ java -jar dist/gtfsplitter.jar -m stack -p 1000 -o TMP ~/jeter.gtf.gz 
INFO	2019-08-07 11:11:09	SortingCollection	Creating merging iterator from 58 files

$ find TMP/ -name "*.gtf"    -exec grep -w -c gene -H '{}' ';'
TMP/af/958f0b968e6c49f3474431cd3aa2df/Group000060.gtf:1000
TMP/ad/7732b2eb185c413c329c1f124bb809/Group000006.gtf:1000
TMP/e0/6a9eb20414e7406ac7a2309ebbe00f/Group000041.gtf:1000
TMP/e0/d2c4d002a926bea4295dbf6db3e831/Group000045.gtf:1000
TMP/f1/d98b194f295a688f1087400ef4cf6b/Group000020.gtf:1000
TMP/f1/07c7ca91843874df6a2c875ea674fd/Group000053.gtf:1000
TMP/6f/ab500a2f6fee4cce61739dd97b11cb/Group000000.gtf:1000
TMP/71/64a6dd2fa9469c1d475d9d6cb78cdc/Group000012.gtf:1000
TMP/a0/3524b44e90561115f6dd0ca3f5a00c/Group000016.gtf:1000
(...)

```



END_DOC
*/
@Program(
		name="gtfsplitter",
		description="Split GTF file per gene, transcript, chromosome...",
		creationDate="20190807",
		modificationDate="20190807",
		keywords= {"gtf"}
		)
public class GtfFileSplitter
	extends Launcher
	{
	private static final Logger LOG = Logger.build(GtfFileSplitter.class).make();
	
	private enum Method {gene,transcript,contig,group /* N-files */,stack /* files contains N genes */};
	
	
	@Parameter(names={"-o","--output"},description= ArchiveFactory.OPT_DESC,required=true)
	private Path outputFile = null;
	@Parameter(names={"-H","--header"},description= "include gtf header")
	private boolean include_header = false;
	@Parameter(names={"-m","--method"},description= "split method. "
			+ "group: will create no more than 'pool-size' gtf files where all the genes will be dispatched. "
			+ "stack: will create files where there won't be more than 'pool-size' gene per file")
	private Method method = Method.gene;
	@Parameter(names={"-p","--pool-size"},description= "size for method=group or method=stack")
	private int pool_size=100;
	@Parameter(names={"-manifest","--manifest"},description="Manifest file containing the path to each gtf")
	private Path manifestFile = null;
	@Parameter(names={"-compress","--compress","--gzip"},description="Gzip output gtf")
	private boolean gzip_gtf =false;

	
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	
		
	public GtfFileSplitter()
		{
		
		}
	
	private static class KeyLine
		{
		final String key;
		final long idx;
		final String line;
		KeyLine(final String key,final long idx,final String line) {
			this.key = key;
			this.idx = idx;
			this.line = line;
			}
		@Override
		public int hashCode() {
			return line.hashCode();
			}
		@Override
		public boolean equals(Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof  KeyLine)) return false;
			return this.line.equals(KeyLine.class.cast(obj).line);
			}
		
		@Override
		public String toString() {
			return this.key;
			}
		}

	private abstract class AbstractSplitter<T extends KeyLine> {
		final GTFCodec codec = new GTFCodec();
		SortingCollection<T> sortingcollection=null;
					
	
		abstract Class<T> getSortedClass();
		abstract SortingCollection.Codec<T> createCodec();
		abstract Comparator<T> createPrimaryComparator();
		abstract Comparator<T> createSecondaryComparator();
		abstract T createKey(String line,long idx,final GTFLine gtfLine);
		abstract String createFilename(T key);
		
		void printGeneLinesIfAny(PrintWriter pw,T key) {
			
		}
		
		int run(final List<String> args) throws IOException {
		BufferedReader in = null;
		CloseableIterator<T> iter=null;
		ArchiveFactory archiveFactory = null;
		boolean got_gtf_record = false;
		final List<String> headerLines = new ArrayList<>();
		PrintWriter manifest = null;
		try {
			
			
			in = GtfFileSplitter.this.openBufferedReader(oneFileOrNull(args));
			archiveFactory = ArchiveFactory.open(GtfFileSplitter.this.outputFile);

			
			this.sortingcollection = SortingCollection.newInstance(
					getSortedClass(),
					createCodec(),
					createPrimaryComparator(),
					GtfFileSplitter.this.writingSortingCollection.getMaxRecordsInRam(),
					GtfFileSplitter.this.writingSortingCollection.getTmpPaths()
					);
			this.sortingcollection.setDestructiveIteration(true);
			long nLines=0;
			for(;;)
				{
				final String line = in.readLine();
				if(line==null) break;
				++nLines;
				if(StringUtils.isBlank(line)) continue;
				if(line.startsWith("#")) {
					if(!got_gtf_record && include_header) {
						headerLines.add(line);
						}
					continue;
					}
				final GTFLine gtfLine =codec.decode(line);
				if(gtfLine==null) continue;
				got_gtf_record = true;
				
				final T key = createKey(line,nLines,gtfLine);
				if(key==null) continue;
				sortingcollection.add(key);
				}
			this.sortingcollection.doneAdding();
			
			manifest = new PrintWriter(GtfFileSplitter.this.manifestFile==null?new NullOuputStream():IOUtils.openPathForWriting(manifestFile));

			Comparator<T> cmp2 = createSecondaryComparator();
			iter = this.sortingcollection.iterator();
			PeekableIterator<T> peekIter = new PeekableIterator<>(iter);
			while(peekIter.hasNext()) {
				final T first = peekIter.next();
				final String basename = createFilename(first);

				
				
				
				final String md5 = StringUtils.md5(basename);
				final String filename =  md5.substring(0,2) + File.separatorChar + md5.substring(2) + File.separator+basename.replaceAll("[/\\:]", "_") + 
						".gtf" + (gzip_gtf?".gz":"");
				
				final PrintWriter os;
				final BlockCompressedOutputStream bcos;
				if(!gzip_gtf) {
					os = archiveFactory.openWriter(filename);
					bcos= null;
					}
				else
					{
					bcos = new BlockCompressedOutputStream
						(archiveFactory.openOuputStream(filename),
						Paths.get(filename)	
						);
					os = new PrintWriter(bcos);
					}
				if(include_header) {
					for(final String h:headerLines) os.println(h);
					}
				printGeneLinesIfAny(os,first);
				os.println(first.line);
				while(peekIter.hasNext()) {
					final T other = peekIter.peek();
					if(cmp2.compare(first,other)!=0) break;
					os.println(peekIter.next().line);
					}
				os.flush();
				os.close();
				
				manifest.println((archiveFactory.isTarOrZipArchive()?"":GtfFileSplitter.this.outputFile.toString()+File.separator)+filename);
				}
			peekIter.close();		
			iter.close();
			manifest.flush();
			manifest.close();
			manifest=null;
			archiveFactory.close();
			return RETURN_OK;
			}
		finally
			{
			CloserUtil.close(manifest);
			if(sortingcollection!=null) sortingcollection.cleanup();
			}
		}
	}
	
	private abstract class SimpleSplitter extends AbstractSplitter<KeyLine> {
		private class MyCodec extends  AbstractDataCodec<KeyLine>
			{
			@Override
			public KeyLine decode(final DataInputStream dis) throws IOException {
				String g;
				try {
					g=dis.readUTF();
				} catch(EOFException err) { return null;}
				final long idx = dis.readLong();
				final String line = dis.readUTF();
				return new KeyLine(g,idx,line);
				}
			@Override
			public void encode(final DataOutputStream dos,final  KeyLine object) throws IOException {
				dos.writeUTF(object.key);
				dos.writeLong(object.idx);
				dos.writeUTF(object.line);
				}
			@Override
			public AbstractDataCodec<KeyLine> clone() {
				return new MyCodec();
				}
			}
		@Override
		String createFilename(KeyLine key) {
			return key.key;
			}
		@Override
		SortingCollection.Codec<KeyLine> createCodec() {
			return new MyCodec();
			}
		@Override
		Comparator<KeyLine> createSecondaryComparator() {
			return (A,B) -> A.key.compareTo(B.key);
			}
		@Override
		Comparator<KeyLine> createPrimaryComparator() {
			final Comparator<KeyLine> cmp=createSecondaryComparator();
			return (A,B) -> {
				int i= cmp.compare(A, B);
				if(i!=0) return i;
				return Long.compare(A.idx, B.idx);
				};
			}
	
		@Override
		Class<KeyLine> getSortedClass() {
			return KeyLine.class;
			}
		}
	
	private  class ContigSplitter extends SimpleSplitter {
		@Override
		KeyLine createKey(String line, long idx,GTFLine gtfLine) {
			return new KeyLine(gtfLine.getContig(),idx,line);
			}
		}
	
	private  class GeneSplitter extends SimpleSplitter {
		@Override
		KeyLine createKey(String line, long idx,GTFLine gtfLine) {
			final String gene_id = gtfLine.getAttribute("gene_id");
			if(StringUtils.isBlank(gene_id)) {
				LOG.warn("Skipping/no gene_id "+line);
				return null;
				}
			
			return new KeyLine(gtfLine.getContig()+"_"+gene_id,idx,line);
			}
		}
	private  class TranscriptSplitter extends SimpleSplitter {
		private final Map<String,String> contig_gene2line = new HashMap<String, String>();
		
		
		
		@Override
		void printGeneLinesIfAny(PrintWriter pw,KeyLine key) {
			GTFLine gtfLine = codec.decode(key.line);
			final String gene_id = gtfLine.getAttribute("gene_id");
			if(StringUtils.isBlank(gene_id)) {
				throw new IllegalStateException();
				}
			String line = contig_gene2line.get(gtfLine.getContig()+"~"+gene_id);
			if(StringUtils.isBlank(line)) {
				LOG.warn("not gene found for "+gtfLine);
				return;
			}
			pw.println(line);
			}
		
		@Override
		KeyLine createKey(String line, long idx,GTFLine gtfLine) {
			final String gene_id = gtfLine.getAttribute("gene_id");
			if(StringUtils.isBlank(gene_id)) {
				LOG.warn("Skipping/no gene_id "+line);
				return null;
				}
			if(gtfLine.getType().equals("gene")) {
				contig_gene2line.put(gtfLine.getContig()+"~"+gene_id,line);
				return null;
				}
			final String transcript_id = gtfLine.getAttribute("transcript_id");
			if(StringUtils.isBlank(transcript_id)) {
				return null;
				}
			
			return new KeyLine(gtfLine.getContig()+"_"+transcript_id,idx,line);
			}
		}
	private  abstract class AbstractGeneGroupSplitter extends SimpleSplitter {
		protected final Map<String,String> gene2pool = new HashMap<>();
		int count_pool = 0;
		
		abstract String getGroupId(String contig,String gene_id);
		
		@Override
		KeyLine createKey(String line, long idx,GTFLine gtfLine) {
			final String gene_id = gtfLine.getAttribute("gene_id");
			if(StringUtils.isBlank(gene_id)) {
				LOG.warn("Skipping/no gene_id "+line);
				return null;
				}
			String group_id = getGroupId(gtfLine.getContig(),gene_id);
			return new KeyLine(group_id,idx,line);
			}
		}
	/** group genes per group into 'pool_size' groups */
	private  class GeneGroupSplitter extends AbstractGeneGroupSplitter {
		@Override
		String getGroupId(String contig, String gene_id) {
			final String key = contig+"~"+gene_id;
			String group = gene2pool.get(key);
			if(group==null) {
				count_pool++;
				group = "Group"+String.format("%06d", count_pool%GtfFileSplitter.this.pool_size);
				gene2pool.put(key,group);
				}
			return group;
			}
		}
	/** group genes per group of 'pool_size' items */
	private  class GeneStackSplitter extends AbstractGeneGroupSplitter {
		int count=0;
		@Override
		String getGroupId(final String contig,final String gene_id) {
			final String key = contig+"~"+gene_id;
			String group = gene2pool.get(key);
			if(group==null) {
			
				if(count +1 > GtfFileSplitter.this.pool_size) {
					count=1;
					count_pool++;
					}
				else
					{
					count++;
					}
				group = "Group"+String.format("%06d",count_pool);
				gene2pool.put(key,group);
				
				}
			return group;
			}
		}
	@Override
	public int doWork(final List<String> args) {
		if(this.pool_size<=0) {
			LOG.error("bad value for pool-size");
			return -1;
		}
		try
			{
			final AbstractSplitter<?> splitter;
			switch(this.method) {
				case contig: splitter = new ContigSplitter();break;
				case gene: splitter = new GeneSplitter();break;
				case transcript: splitter = new TranscriptSplitter();break;
				case stack: splitter = new GeneStackSplitter();break;
				case group: splitter = new GeneGroupSplitter();break;
				default: throw new IllegalStateException("method: "+this.method);
				}
			
			return splitter.run(args);
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	 	
	
	public static void main(final String[] args)
		{
		new GtfFileSplitter().instanceMainWithExit(args);
		}
	}
