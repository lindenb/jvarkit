/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.bedclustername;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLine;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.FileExtensions;

/**
BEGIN_DOC

## Example

```
$ gunzip -c src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz | awk -F '\t' '$3=="exon"' | \
	java -jar dist/jvarkit.jar gtf2bed -c 'gene_name' |\
	java -jar dist/jvarkit.jar bedclustername -o jeter.zip  --size 100 --merge 1

$ unzip -l jeter.zip 
Archive:  jeter.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
     7332  2025-04-28 14:37   cluster.000000001.bed
     1740  2025-04-28 14:37   cluster.000000002.bed
     1456  2025-04-28 14:37   cluster.000000003.bed
---------                     -------
    10528                     3 files


$ unzip -p jeter.zip cluster.000000002.bed | head
22	41697718	41697776	ZC3H7B
22	41716664	41716717	ZC3H7B
22	41721567	41721601	ZC3H7B
22	41721724	41721922	ZC3H7B
22	41723209	41723368	ZC3H7B
22	41726026	41726107	ZC3H7B
22	41728174	41732847	ZC3H7B
22	41734316	41734359	ZC3H7B
22	41735004	41735195	ZC3H7B
22	41735819	41736141	ZC3H7B

```


END_DOC

 */
@Program(
		name="bedclustername",
		description="Clusters a BED file into a set of BED files using the 4th column of the bed name.",
		keywords={"bed","chromosome","contig"},
		creationDate="2050428",
		modificationDate="2050428",
		jvarkit_amalgamion =  true,
		menu="BED Manipulation"
		)
public class BedClusterName
	extends Launcher
	{
	private static final Logger LOG = Logger.of(BedClusterName.class);
	private enum FormatOut {BED,BED_GZ,INTERVAL_LIST,INTERVAL_LIST_GZ}
	@Parameter(names={"-o","--out"},description=ArchiveFactory.OPT_DESC,required=true)
	private Path outputFile= null;
	@Parameter(names={"-F","--format"},description=ArchiveFactory.OPT_DESC)
	private FormatOut outputFormat = FormatOut.BED;
	@Parameter(names={"-J","--jobs"},description="number of clusters. (or specify --size or --window-size/--window-shif)")
	private int number_of_jobs=-1;
	@Parameter(names={"-S","--size"},description="number of bases max per bin. (or specify --jobs or --window-size/--window-shif). "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.LongStringConverter.class,splitter=NoSplitter.class)
	private long long_length_per_bin=-1L;
	@Parameter(names={"-C","--contig","--chromosome"},description="group by chromosome.")
	private boolean group_by_contig = false;
	@Parameter(names={"-R","--reference"},description="For sorting and /or writing interval_list," +DICTIONARY_SOURCE)
	private Path refFaidx = null;
	@Parameter(names={"-m","--manifest"},description="Manifest Bed file output containing chrom/start/end of each gene")
	private Path manifestFile = null;
	@Parameter(names={"--md5-dir","--sub-dir"},description="prevent the creation of too many files in the same directory. Create some intermediate directories based on filename's md5.")
	private boolean create_md5_dir_flag=false;
	@Parameter(names={"--merge-distance","--merge"},description="if greater than -1 Merge overlapping record for the same name within a distance of 'x'. " + DistanceParser.OPT_DESCRIPTION, converter = DistanceParser.StringConverter.class)
	private int merge_distance=-1;
	
	private int id_generator =0;
	
	private static class Batch {
		final String name;
		private final List<BedLine> beds = new ArrayList<>();
		Batch(final String name) {
			this.name=name;
			}
		long getLengthOnReference() {
			return beds.stream().mapToLong(B->B.getLengthOnReference()).sum();
		}

		void merge(final int distance) {
			Collections.sort(beds,(A,B)->Integer.compare(A.getStart(), B.getEnd()));
			int i=0;
			while(i+1 < beds.size()) {
				final BedLine bed0= beds.get(i+0);
				final BedLine bed1= beds.get(i+1);
				if(bed0.withinDistanceOf(bed1, distance)) {
					final String[] tokens = bed0.toStringArray();
					tokens[1]=String.valueOf(Math.min(bed0.getStart(), bed1.getStart()));
					tokens[2]=String.valueOf(Math.max(bed0.getEnd(), bed1.getEnd()));
					beds.set(i,new BedLine(tokens));
					beds.remove(i+1);
					}
				else
					{
					i++;
					}
				}
			}

		public String getUniqContig() {
			final Set<String> contigs = this.beds.stream().map(B->B.getContig()).collect(Collectors.toSet());
			if(contigs.size()!=1) {
				throw new IllegalArgumentException("Cannot group by contig when there is more than one contig per batch "+this.beds);
				}
			return contigs.iterator().next();
			}

		@Override
		public String toString() {
			return beds.toString();
			}
	}

	private static class Cluster extends AbstractList<Batch>{
		private long sum_length=0L;
		private final List<Batch> intervals = new ArrayList<>();
		
		@Override
		public boolean add(final Batch si) {
			this.intervals.add(si);
			this.sum_length+=si.getLengthOnReference();
			return true;
		}
		
		/** get sum of lengths */
		long getSumLength(final Batch malus) {
			return this.sum_length + (malus==null?0:malus.getLengthOnReference());
		}
		/** get average length in cluster */
		double getAvgLength() {
			final int n= this.intervals.size();
			final double total= getSumLength(null);
			return total/n;
		}
		
		/** get standard deviation of intervals length if 'malus' is added to the pool */
		double getStdDev() {
			final double avg = this.getAvgLength();
			final int n=intervals.size();
			
			return  (
					intervals.stream().mapToDouble(L->Math.abs(avg-L.getLengthOnReference())).sum()
					) / n;
			}
		
		@Override
		public Batch get(int index) {
			return this.intervals.get(index);
			}
		
		@Override
		public int size() {
			return this.intervals.size();
		}
	}
	
	
	
	private final Comparator<BedLine> defaultIntervalCmp=(B1,B2)->{
			int i = B1.getContig().compareTo(B2.getContig());
			if(i!=0) return i;
			i = Integer.compare(B1.getStart(), B2.getStart());
			if(i!=0) return i;
			i = Integer.compare(B1.getEnd(), B2.getEnd());
			return i;
			};

	private String md5dir(final String s) {
		if(!this.create_md5_dir_flag) return s;
		final String md5 = StringUtils.md5(s);
		return String.join(File.separator, md5.substring(0,2),md5.substring(2,4),s);
		}
			
	
	private long apply_cluster(
			final PrintWriter manifest,
			final ArchiveFactory archiveFactory,
			final List<Batch> src,
			final Comparator<BedLine> sorter,
			final SAMSequenceDictionary dict
			) throws IOException
		{
		long length_out=0L;
		final List<Cluster> clusters = new ArrayList<>(Math.max(this.number_of_jobs,100));
		final LinkedList<Batch> list = new LinkedList<>(src);
		
	
		//sort by decreasing size
		Collections.sort(list,(B1,B2)->Long.compare(B2.getLengthOnReference(),B1.getLengthOnReference()));
			
		
		/** N-CLUSTER ALGORITHM *************************************/
		if(this.number_of_jobs>0) {
			while(!list.isEmpty()) {
				final Batch first=list.pop();
				if(clusters.size() < this.number_of_jobs) {
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
		/** LENGTH ALGORITHM *************************************/
		else // group by size
			{
			while(!list.isEmpty()) {
				final Batch first=list.pop();
				int y=0;
				while(y<clusters.size()) {
					final Cluster cluster = clusters.get(y);
					if(cluster.getSumLength(first) <= this.long_length_per_bin) {
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
				
		for(final Cluster cluster : clusters) {
			final List<BedLine> all_beds = cluster.intervals.stream().
					flatMap(C->C.beds.stream()).
					sorted(sorter).
					collect(Collectors.toList());
			final String prefix;
			if(this.group_by_contig) {
				final int chromStart = all_beds.stream().mapToInt(R->R.getStart()-1).min().orElse(0);
				final int chromEnd = all_beds.stream().mapToInt(R->R.getEnd()).max().orElse(0);
				manifest.print(all_beds.get(0).getContig());
				manifest.print("\t");
				manifest.print(chromStart);
				manifest.print("\t");
				manifest.print(chromEnd);
				prefix= all_beds.get(0).getContig()+"_"+(chromStart+1)+"_"+chromEnd;
				}
			else if(cluster.intervals.size()==1) {
				prefix="cluster."+cluster.intervals.get(0).name.replaceAll("[^A-Z_a-z0-9]","_");
				}
			else
				{
				prefix="cluster";
				}

			boolean save_as_interval_list;
			final String suffix;
			switch(this.outputFormat) {
				case INTERVAL_LIST:
				case INTERVAL_LIST_GZ:
					{
					suffix = FileExtensions.INTERVAL_LIST  +(this.outputFormat.name().endsWith("GZ")?".gz":"");
					save_as_interval_list=true;
					break;
					}
				case BED:
				case BED_GZ:
					{
					suffix = FileExtensions.BED +(this.outputFormat.name().endsWith("GZ")?".gz":"");
					save_as_interval_list=false;
					break;
					}
				default: throw new IllegalArgumentException("err");
				}
			
			final String filename = md5dir(String.format("%s.%09d%s",prefix, ++id_generator,suffix));

			manifest.print(archiveFactory.isTarOrZipArchive()?filename:this.outputFile.resolve(filename).toAbsolutePath());
			manifest.print("\t");
			manifest.print(cluster.size());
			manifest.print("\t");
			manifest.print(cluster.getSumLength(null));
			manifest.print("\t");
			manifest.print((int)cluster.getAvgLength());
			manifest.print("\t");
			manifest.print((int)cluster.getStdDev());
			manifest.print("\t");
			manifest.print(String.join(",", cluster.intervals.stream().map(B->B.name).collect(Collectors.toCollection(TreeSet::new))));
			manifest.println();

			try(final OutputStream ps = archiveFactory.openOuputStream(filename)) {
				final PrintWriter pw;
				final BlockCompressedOutputStream bcos;
				if(this.outputFormat.name().endsWith("GZ")) {
					bcos = new BlockCompressedOutputStream(ps, (Path)null);
					pw=new PrintWriter(bcos);
				} else
					{
					bcos = null;
					pw=new PrintWriter(ps);
					}
				
				if(save_as_interval_list) {
					final SAMFileHeader samHeader = new SAMFileHeader(dict);
					samHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
					final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
					JVarkitVersion.getInstance().addMetaData(this, samHeader);
					codec.encode(pw, samHeader);
					}
				
					
					
				for(final BedLine r: all_beds) {
					pw.print(r.getContig());
					pw.print("\t");
					pw.print(r.getStart()-(save_as_interval_list?0:1));//convert to bed if needed
					pw.print("\t");
					pw.print(r.getEnd());
					for(int c=3;c< r.getColumnCount();c++) {
						pw.print("\t");
						pw.print(r.getOrDefault(c, "."));
						}
					pw.println();
					length_out+= r.getLengthOnReference();
					}
					
				pw.flush();
				pw.close();
				if(bcos!=null) {bcos.flush();bcos.close();}
				ps.flush();
				}
			}
		return length_out;
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.outputFormat.name().startsWith("INTERVAL") && this.refFaidx==null) {
			LOG.error("REF must be specified when saving as interval list");
			return -1;
			}
		
		
		final int n_args = (this.number_of_jobs>0?1:0)+
				(this.long_length_per_bin>0?1:0)
				;
		
		if(n_args==0) {
			LOG.error("at least --jobs or --size  must be specified.");
			return -1;
			}
		if(n_args>1) {
			LOG.error(" --jobs OR --size must be specified. Not both.");
			return -1;
			}

		
		
		try
			{
			final SAMSequenceDictionary dict;
			final ContigNameConverter converter;
			final Comparator<String> contigComparator;
			final Comparator<BedLine> intervalComparator;
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
			
			
			
			try(ArchiveFactory archiveFactory = ArchiveFactory.open(this.outputFile)) {
				if(this.outputFormat.name().endsWith("GZ") && archiveFactory.isTarOrZipArchive()) {
					archiveFactory.setCompressionLevel(0);
					}
				
				try(PrintWriter manifest = new PrintWriter(this.manifestFile==null?new NullOuputStream():IOUtils.openPathForWriting(this.manifestFile))) {
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
					manifest.print("\t");
					manifest.print("names");
					manifest.println();

					final String input = oneFileOrNull(args);
					try(BedLineReader br = new BedLineReader(super.openBufferedReader(input),input)) {
						br.setValidationStringency(ValidationStringency.LENIENT);
						br.setContigNameConverter(converter);
						final Map<String,Batch> name2batch = new HashMap<>();
						
						while(br.hasNext()) {
							final BedLine bed = br.next();
							final String name = bed.getOrDefault(3, "");
							if(StringUtils.isBlank(name)) {
								LOG.error("empty 4th column in "+bed);
								return -1;
								}
							Batch batch  = name2batch.get(name);
							if(batch==null) {
								batch=new Batch(name);
								name2batch.put(name, batch);
								}
							batch.beds.add(bed);
							}

						if(this.merge_distance>=0) {
							for(Batch batch:name2batch.values()) {
								batch.merge(this.merge_distance);
								}
							}

						List<List<Batch>> list_of_batches = Collections.singletonList(new ArrayList<>(name2batch.values()));
						
						
						final long length_in = name2batch.values().stream().flatMap(it->it.beds.stream()).mapToLong(L->L.getLengthOnReference()).sum();
						
						if(this.group_by_contig) {
							final List<List<Batch>> list2= new ArrayList<>();
							for(List<Batch> l1:list_of_batches) {
								final Map<String,List<Batch>> contig2lines = new TreeMap<>(contigComparator);
								contig2lines.putAll(l1.stream().collect(Collectors.groupingBy(B->B.getUniqContig())));
								for(final String ctg:contig2lines.keySet()) {
									list2.add(new ArrayList<>(contig2lines.get(ctg)));
									}
								}
							list_of_batches = list2;
							}
						
						long length_out = 0L;
						for(final List<Batch> l1:list_of_batches) {
							length_out += apply_cluster(manifest,archiveFactory,l1,intervalComparator,dict);
							}
						
						// let me be paranoid
						if(length_in!=length_out) {
							throw new IllegalStateException("length-in!=length-out");
							}
						}
					manifest.flush();
					}
				}
			
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		}
	

	public static void main(final String[] args)
		{
		new BedClusterName().instanceMainWithExit(args);
		}
	}
