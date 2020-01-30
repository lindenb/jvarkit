/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloserUtil;

/**
BEGIN_DOC



END_DOC

 */
@Program(
		name="bedcluster",
		description="Convert the names of the chromosomes in a Bed file",
		keywords={"bed","chromosome","contig"},
		creationDate="20200130",
		modificationDate="20200130"
		)
public class BedCluster
	extends Launcher
	{
	private static final Logger LOG = Logger.build(BedCluster.class).make();
	
	@Parameter(names={"-o","--out"},description=ArchiveFactory.OPT_DESC,required=true)
	private Path outputFile= null;
	@Parameter(names={"-J","--jobs"},description="number of clusters.")
	private int njobs=100;
	@Parameter(names={"-C","--cgonti"},description="group by chromosome.")
	private boolean group_by_contig = false;
	@Parameter(names={"-R","--reference"},description="For sorting. "+INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path refFaidx = null;
	@Parameter(names={"--merge"},description="Merge bed records")
	private boolean merge_bed_records=false;
	@Parameter(names={"-m","--manifest"},description="Manifest Bed file output containing chrom/start/end of each gene")
	private Path manifestFile = null;
	@Parameter(names={"--compress"},description="Compress with bgzi")
	private boolean do_compress=false;

	private int id_generator =0;
	
	private double getStdDev(final List<SimpleInterval> list,final SimpleInterval malus) {
		final int n=list.size()+1;
		final double total= list.stream().mapToDouble(L->L.getLengthOnReference()).sum() + malus.getLengthOnReference();
		final double avg = total/n;
		
		return  (Math.abs(avg-malus.getLengthOnReference()) + list.stream().mapToDouble(L->Math.abs(avg-L.getLengthOnReference())).sum()) / n;
		
	}
	
	private final Comparator<SimpleInterval> defaultIntervalCmp=(B1,B2)->{
			int i = B1.getContig().compareTo(B2.getContig());
			if(i!=0) return i;
			i = Integer.compare(B1.getStart(), B2.getStart());
			if(i!=0) return i;
			i = Integer.compare(B1.getEnd(), B2.getEnd());
			return i;
			};
	
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
	
	private void apply_cluster(final PrintWriter manifest,final ArchiveFactory archiveFactory,final List<SimpleInterval> src)
			throws IOException
		{
		final List<List<SimpleInterval>> clusters = new ArrayList<>(this.njobs);
		final LinkedList<SimpleInterval> list = new LinkedList<>(mergeBedRecords(src));
		//sort by decreasing size
		Collections.sort(list,(B1,B2)->Integer.compare(B2.getLengthOnReference(),B1.getLengthOnReference()));
		
		while(!list.isEmpty()) {
			final SimpleInterval first=list.pop();
			if(clusters.size()<this.njobs) {
				final List<SimpleInterval> c = new ArrayList<>();
				c.add(first);
				clusters.add(c);
			} else
				{
				int best_idx=-1;
				double best_stddev=-1;
				for(int y=0;y< clusters.size();++y) {
					
					final double stdev= getStdDev(clusters.get(y),first);
					if(best_idx==-1 ||stdev<best_stddev ) {
						best_idx=y;
						best_stddev = stdev;
						}
					}
				clusters.get(best_idx).add(first);
				}
			}
		
		
		
		for(final List<SimpleInterval> cluster : clusters) {
			Collections.sort(cluster,defaultIntervalCmp);
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
			final String filename = String.format(prefix + ".%05d.bed%s", ++id_generator,do_compress?".gz":"");

			manifest.print(filename);
			manifest.print("\t");
			manifest.print(cluster.size());
			manifest.println();
			
			
			final OutputStream ps = archiveFactory.openOuputStream(filename);
			final PrintWriter pw;
			final BlockCompressedOutputStream bcos;
			if(!this.do_compress) {
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
	
	@Override
	public int doWork(final List<String> args) {
		final BedLineCodec codec = new BedLineCodec();
		ArchiveFactory archiveFactory=null;
		PrintWriter manifest=null;

		try
			{
			final ContigNameConverter converter;
			final Comparator<String> contigComparator;
			if(this.refFaidx==null) {
				converter = ContigNameConverter.getIdentity();
				contigComparator = (A,B)->A.compareTo(B);
			} else
			{
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.refFaidx);
				converter = ContigNameConverter.fromOneDictionary(dict);
				contigComparator = new ContigDictComparator(dict);
			}
			
			archiveFactory = ArchiveFactory.open(this.outputFile);
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
						apply_cluster(manifest,archiveFactory,contig2lines.get(ctg));
						}
					}
				else
					{
					apply_cluster(manifest,archiveFactory,st1.collect(Collectors.toList()));
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
