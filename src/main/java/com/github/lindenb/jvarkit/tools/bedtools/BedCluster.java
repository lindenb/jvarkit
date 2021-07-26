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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.bedtools;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;

/**
BEGIN_DOC

## Example

```
$ java -jar dist/bedcluster.jar -j 10 -m jeter.mf -o jeter.zip --compress --contig test.bed

$ head jeter.mf

#chrom	start	end	filename	number_of_records	sum-length	avg-length	stddev-size
1	57460450	59012406	1_57460451_59012406.000000001.bed.gz	1	1551956	1551956	0
1	48998526	50489585	1_48998527_50489585.000000002.bed.gz	1	1491059	1491059	0
1	901876	248790491	1_901877_248790491.000000003.bed.gz	49	1393594	28440	41317
1	1567473	248154506	1_1567474_248154506.000000004.bed.gz	50	1393602	27872	39042
1	470970	229841608	1_470971_229841608.000000005.bed.gz	51	1393594	27325	37502
1	160445	248041507	1_160446_248041507.000000006.bed.gz	51	1393602	27325	37487
1	5647427	246685894	1_5647428_246685894.000000007.bed.gz	51	1393601	27325	37433
1	34553	245778447	1_34554_245778447.000000008.bed.gz	51	1393601	27325	37186
1	696290	247495148	1_696291_247495148.000000009.bed.gz	52	1393594	26799	35931

$ unzip -l jeter.zip | head
Archive:  jeter.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
       76  2020-01-31 12:45   1_57460451_59012406.000000001.bed.gz
       74  2020-01-31 12:45   1_48998527_50489585.000000002.bed.gz
      487  2020-01-31 12:45   1_901877_248790491.000000003.bed.gz
      494  2020-01-31 12:45   1_1567474_248154506.000000004.bed.gz
      494  2020-01-31 12:45   1_470971_229841608.000000005.bed.gz
      511  2020-01-31 12:45   1_160446_248041507.000000006.bed.gz
      511  2020-01-31 12:45   1_5647428_246685894.000000007.bed.gz

 unzip -p jeter.zip 1_5647428_246685894.000000007.bed.gz | gunzip  -c  |head
1	5647427	5728355
1	8921060	8939308
1	11128527	11133154
1	13216355	13219581
1	16787442	16794976
1	20465818	20476879
1	21737952	21739786
1	26145130	26159432
1	37627163	37627235
1	38021842	38022108
```


END_DOC

 */
@Program(
		name="bedcluster",
		description="Clusters a BED file into a set of BED files.",
		keywords={"bed","chromosome","contig"},
		creationDate="20200130",
		modificationDate="20210726",
		biostars= {424828}
		)
public class BedCluster
	extends Launcher
	{
	private static final Logger LOG = Logger.build(BedCluster.class).make();
	
	@Parameter(names={"-o","--out"},description=ArchiveFactory.OPT_DESC,required=true)
	private Path outputFile= null;
	@Parameter(names={"-J","--jobs"},description="number of clusters. (or specify --size)")
	private int number_of_jobs=-1;
	@Parameter(names={"-S","--size"},description="number of bases max per bin. (or specify --jobs). "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.LongStringConverter.class,splitter=NoSplitter.class)
	private long long_length_per_bin=-1L;
	@Parameter(names={"-C","--contig","--chromosome"},description="group by chromosome.")
	private boolean group_by_contig = false;
	@Parameter(names={"-R","--reference"},description="For Sorting." +DICTIONARY_SOURCE)
	private Path refFaidx = null;
	@Parameter(names={"--merge"},description="Merge overlapping bed records before clustering")
	private boolean merge_bed_records=false;
	@Parameter(names={"-m","--manifest"},description="Manifest Bed file output containing chrom/start/end of each gene")
	private Path manifestFile = null;
	@Parameter(names={"--compress"},description="Compress bed with bgz.")
	private boolean do_compress=false;
	@Parameter(names={"--interval-list"},description="Save as htsjdk interval list (https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/util/IntervalList.html) instead of bed.")
	private boolean save_as_interval_list=false;
	@Parameter(names={"--consecutive"},description="When using option --size only use consecutive ordered region. Default is to find the best region anywhere.")
	private boolean consecutive_flags = false;

	
	
	private int id_generator =0;
	
	private static class Cluster extends AbstractList<SimpleInterval>{
		private long sum_length=0L;
		private final List<SimpleInterval> intervals = new ArrayList<>();
		
		@Override
		public boolean add(final SimpleInterval si) {
			this.intervals.add(si);
			this.sum_length+=si.getLengthOnReference();
			return true;
		}
		
		/** get sum of lengths */
		long getSumLength(final SimpleInterval malus) {
			return this.sum_length + (malus==null?0:malus.getLengthOnReference());
		}
		/** get average length in cluster */
		double getAvgLength(final SimpleInterval malus) {
			final int n= this.intervals.size()+(malus==null?0:1);
			final double total= getSumLength(malus);
			return total/n;
		}
		
		/** get standard deviation of intervals length if 'malus' is added to the pool */
		double getStdDev(final SimpleInterval malus) {
			final double avg = this.getAvgLength(malus);
			final int n=intervals.size()+(malus==null?0:1);
			
			return  (
					(malus==null?0.0:Math.abs(avg-malus.getLengthOnReference())) + 
					intervals.stream().mapToDouble(L->Math.abs(avg-L.getLengthOnReference())).sum()
					) / n;
			}
		
		@Override
		public SimpleInterval get(int index) {
			return this.intervals.get(index);
			}
		
		@Override
		public int size() {
			return this.intervals.size();
		}
	}
	

	
	

	
	
	
	private final Comparator<SimpleInterval> defaultIntervalCmp=(B1,B2)->{
			int i = B1.getContig().compareTo(B2.getContig());
			if(i!=0) return i;
			i = Integer.compare(B1.getStart(), B2.getStart());
			if(i!=0) return i;
			i = Integer.compare(B1.getEnd(), B2.getEnd());
			return i;
			};
	
	/** merge overlapping bed records */
	private List<SimpleInterval> mergeBedRecords(final List<SimpleInterval> src) {
		if(!this.merge_bed_records) return src;
		final List<SimpleInterval> list = new ArrayList<>(src);
		Collections.sort(list,defaultIntervalCmp);
		int i=0;
		while(i +1< list.size()) {
			final SimpleInterval r1 = list.get(i);
			final SimpleInterval r2 = list.get(i+1);
			if(r1.overlaps(r2)) {
				list.remove(i+1);
				list.set(i, r1.merge(r2));
			} else
			{
				i++;
			}
		}
	return list;
	}
	
	private void apply_cluster(
			final PrintWriter manifest,
			final ArchiveFactory archiveFactory,
			final List<SimpleInterval> src,
			final Comparator<SimpleInterval> sorter,
			final SAMSequenceDictionary dict)
			throws IOException
		{
		final List<Cluster> clusters = new ArrayList<>(Math.max(this.number_of_jobs,100));
		final LinkedList<SimpleInterval> list = new LinkedList<>(mergeBedRecords(src));
		//sort by decreasing size
		Collections.sort(list,(B1,B2)->Integer.compare(B2.getLengthOnReference(),B1.getLengthOnReference()));
		
		if(number_of_jobs>0) {
			while(!list.isEmpty()) {
				final SimpleInterval first=list.pop();
				if(clusters.size()<this.number_of_jobs) {
					final Cluster c = new Cluster();
					c.add(first);
					clusters.add(c);
					}
				else
					{
					int best_idx=-1;
					double best_length=-1;
					for(int y=0;y< clusters.size();++y) {
						final double total_length = clusters.get(y).getSumLength(first);
						if(best_idx==-1 ||total_length<best_length ) {
							best_idx=y;
							best_length = total_length;
							}
						}
					clusters.get(best_idx).add(first);
					}
				}
			}
		else if( this.consecutive_flags) {// group by size using consecutive regions
			while(!list.isEmpty()) {
				final Cluster cluster = new Cluster();
				cluster.add(list.pop());
				while(!list.isEmpty()) {
					final SimpleInterval si = list.peek();
					if(cluster.getSumLength(si) > this.long_length_per_bin) break;
					cluster.add(list.pop());
					}
				clusters.add(cluster);
				}
			}
		else // group by size
			{
			while(!list.isEmpty()) {
				final SimpleInterval first=list.pop();
				int y=0;
				while(y<clusters.size()) {
					final Cluster cluster = clusters.get(y);
					if(cluster.getSumLength(first)<=this.long_length_per_bin) {
						cluster.add(first);
						break;
						}
					y++;
					}
				if(y==clusters.size()) {
					final Cluster cluster = new Cluster();
					cluster.add(first);
					clusters.add(cluster);
					}
				}
			}
		
		final Path tmpPath = (this.save_as_interval_list?Files.createTempFile("data.", FileExtensions.INTERVAL_LIST):null);
		
		for(final Cluster cluster : clusters) {
			Collections.sort(cluster.intervals,sorter);
			final String prefix;

			if(this.group_by_contig) {
				final int chromStart = cluster.stream().mapToInt(R->R.getStart()-1).min().orElse(0);
				final int chromEnd = cluster.stream().mapToInt(R->R.getEnd()).max().orElse(0);
				manifest.print(cluster.get(0).getContig());
				manifest.print("\t");
				manifest.print(chromStart);
				manifest.print("\t");
				manifest.print(chromEnd);
				manifest.print("\t");
				prefix= cluster.get(0).getContig()+"_"+(chromStart+1)+"_"+chromEnd;
			} else
				{
				prefix="cluster";
				}
			
			final String suffix;
			if(this.save_as_interval_list) {
				suffix = FileExtensions.INTERVAL_LIST;
			} else
				{
				suffix = FileExtensions.BED +(do_compress?".gz":"");
				}
			
			final String filename = String.format("%s.%09d%s",prefix, ++id_generator,suffix);

			manifest.print(archiveFactory.isTarOrZipArchive()?filename:this.outputFile.resolve(filename).toAbsolutePath());
			manifest.print("\t");
			manifest.print(cluster.size());
			manifest.print("\t");
			manifest.print(cluster.getSumLength(null));
			manifest.print("\t");
			manifest.print((int)cluster.getAvgLength(null));
			manifest.print("\t");
			manifest.print((int)cluster.getStdDev(null));

			manifest.println();
			
			if(this.save_as_interval_list) {
				final IntervalList iL = new IntervalList(Objects.requireNonNull(dict));
				iL.addall(cluster.stream().map(R->new Interval(R)).collect(Collectors.toList()));
				iL.write(tmpPath);
				archiveFactory.copyTo(tmpPath, filename);
				Files.deleteIfExists(tmpPath);
				}
			else
				{
				final OutputStream ps = archiveFactory.openOuputStream(filename);
				final PrintWriter pw;
				final BlockCompressedOutputStream bcos;
				if(this.do_compress) {
					bcos = new BlockCompressedOutputStream(ps, (Path)null);
					pw=new PrintWriter(bcos);
				} else
					{
					bcos = null;
					pw=new PrintWriter(ps);
					}
				
				for(final SimpleInterval r: cluster) {
					pw.print(r.getContig());
					pw.print("\t");
					pw.print(r.getStart()-1);
					pw.print("\t");
					pw.print(r.getEnd());
					pw.println();
					}
				pw.flush();
				pw.close();
				if(bcos!=null) {bcos.flush();bcos.close();}
				ps.close();
				}
			}
		
		}
	
	@Override
	public int doWork(final List<String> args) {
		final BedLineCodec codec = new BedLineCodec();
		ArchiveFactory archiveFactory=null;
		PrintWriter manifest=null;
		if(save_as_interval_list && this.refFaidx==null) {
			LOG.error("REF must be specified when saving as interval list");
			return -1;
			}
		if(this.number_of_jobs>0 && this.consecutive_flags) {
			LOG.error("--consecutive cannot be set when --jobs is used.");
			return -1;
			}
		if(this.number_of_jobs<1 && this.long_length_per_bin<1L) {
			LOG.error("at least --jobs or --size must be specified.");
			return -1;
			}
		if(this.number_of_jobs>0 &&  this.long_length_per_bin>0) {
			LOG.error(" --jobs OR --size must be specified. Not both.");
			return -1;
			}

		try
			{
			final SAMSequenceDictionary dict;
			final ContigNameConverter converter;
			final Comparator<String> contigComparator;
			final Comparator<SimpleInterval> intervalComparator;
			if(this.refFaidx==null) {
				dict = null;
				converter = ContigNameConverter.getIdentity();
				contigComparator = (A,B)->A.compareTo(B);
				intervalComparator = this.defaultIntervalCmp;
			} else
			{
				dict = SequenceDictionaryUtils.extractRequired(this.refFaidx);
				converter = ContigNameConverter.fromOneDictionary(dict);
				final ContigDictComparator cmp2 =  new ContigDictComparator(dict);
				contigComparator = cmp2;
				intervalComparator = cmp2.createLocatableComparator();
			}
			
			archiveFactory = ArchiveFactory.open(this.outputFile);
			if(!this.save_as_interval_list && this.do_compress && archiveFactory.isTarOrZipArchive()) {
				archiveFactory.setCompressionLevel(0);
			}
			
			manifest = new PrintWriter(this.manifestFile==null?new NullOuputStream():IOUtils.openPathForWriting(this.manifestFile));
			manifest.print("#");
			if(this.group_by_contig) {
				manifest.print("chrom");
				manifest.print("\t");
				manifest.print("start");
				manifest.print("\t");
				manifest.print("end");
				manifest.print("\t");
			}
			manifest.print("filename");
			manifest.print("\t");
			manifest.print("number_of_records");
			manifest.print("\t");
			manifest.print("sum-length");
			manifest.print("\t");
			manifest.print("avg-length");
			manifest.print("\t");
			manifest.print("stddev-size");
			manifest.println();
			
			try(BufferedReader br = super.openBufferedReader(oneFileOrNull(args))) {
				final Stream<SimpleInterval> st1 = br.lines().
					filter(L->!BedLine.isBedHeader(L)).
					map(L->codec.decode(L)).
					filter(B->B!=null).
					filter(B->!StringUtils.isBlank(converter.apply(B.getContig()))).
					map(B->new SimpleInterval(converter.apply(B.getContig()), B.getStart(), B.getEnd()));
					;
				
				if(this.group_by_contig) {
					final Map<String,List<SimpleInterval>> contig2lines = new TreeMap<>(contigComparator);
					contig2lines.putAll(st1.collect(Collectors.groupingBy(B->B.getContig())));
					for(final String ctg:contig2lines.keySet()) {
						apply_cluster(manifest,archiveFactory,contig2lines.get(ctg),intervalComparator,dict);
						}
					}
				else
					{
					apply_cluster(manifest,archiveFactory,st1.collect(Collectors.toList()),intervalComparator,dict);
					}
				}
			manifest.flush();
			manifest.close();
			archiveFactory.close();
			
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(manifest);
			CloserUtil.close(archiveFactory);
			}
		}
	

	public static void main(final String[] args)
		{
		new BedCluster().instanceMainWithExit(args);
		}
	}
