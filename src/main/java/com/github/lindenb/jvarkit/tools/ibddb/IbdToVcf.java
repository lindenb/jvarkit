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
package com.github.lindenb.jvarkit.tools.ibddb;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;

/**
BEGIN_DOC


````
$ java -jar dist/ibd2vcf.jar  --fam TSF.fam --bed TSF.bed -R human_b37.dict  TSF_Chr22.ibdtxt  |\
	bcftool view -O b -o out.bcf
```

END_DOC
*/
@Program(name="ibd2vcf",
	description="IBD data to VCF",
	keywords={"ibd","vcf"},
	creationDate="20210701",
	modificationDate="20210706",
	generate_doc=false
	)
public class IbdToVcf extends Launcher {
	private static final Logger LOG = Logger.build(IbdToVcf.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference","--dict"},description= DICTIONARY_SOURCE,required=true)
	private Path dictPath = null;
	@Parameter(names={"-B","--bed","--markers"},description="Marker file chrom(tab)start(tab)end(tab)name",required=true)
	private Path markerBedPath = null;
	@Parameter(names={"-P","--fam","--pedigree"},description="Pedigree file. family-id(space)sample-id.",required=true)
	private Path famPath = null;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();

	private SAMSequenceDictionary dictionary;
	private final List<Marker> all_markers = new ArrayList<>();
	private final List<Pair> all_pairs = new ArrayList<>();
	private final Map<Pair,Pair> pair2pair = new HashMap<>();
	private final List<Sample> all_samples = new ArrayList<>();
	private final Pattern SPACES = Pattern.compile("[\\s]+");
	
	private class Sample implements Comparable<Sample> {
		final String family;
		final String name;
		Sample(final String family,final String name) {
			this.family = family;
			this.name = name;
		}
		@Override
		public int compareTo(Sample o) {
			int i= this.family.compareTo(o.family);
			if(i!=0) return i;
			return name.compareTo(o.name);
			}
		@Override
		public String toString() {
			return family+"_"+name;
			}
		}
	
	private class Pair {
		final Sample sample1;
		final Sample sample2;
		int index = -1;
		Pair(final Sample sn1,final Sample sn2) {
			if(sn1.compareTo(sn2)<=0) {
				this.sample1 = sn1;
				this.sample2 = sn2;
				} else {
				this.sample1 = sn2;
				this.sample2 = sn1;
				}
			}
		@Override
		public boolean equals(Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof Pair)) return false;
			final Pair other = Pair.class.cast(obj);
			return this.sample1.equals(other.sample1) &&
					this.sample2.equals(other.sample2)
					;
			}
		@Override
		public int hashCode() {
			return this.sample1.hashCode()*31 + this.sample2.hashCode();
			}
		String getName() {
			return sample1+"|"+sample2;
			}
		@Override
		public String toString() {
			return getName();
			}
	}
	
	
	private static class Data {
		int pair_id;
		int marker_id;
		float[] ibd = new float[3];
		
		public int compare2(final Data d) {
			return Integer.compare(this.marker_id, d.marker_id);
			}
		public int compare1(final Data d) {
			int i= compare2(d);
			if(i!=0) return i;
			return Integer.compare(this.pair_id, d.pair_id);
			}
		}
	
	
	private static class DataCodec extends AbstractDataCodec<Data> {
		@Override
		public Data decode(DataInputStream dis) throws IOException {
			final Data d = new Data();
			try {
				d.pair_id = dis.readInt();
			} catch(EOFException err) {
				return null;
				}
			d.marker_id = dis.readInt();
			for(int i=0;i<d.ibd.length;i++) {
				d.ibd[i]= dis.readFloat();
				}
			return d;
			}
		@Override
		public void encode(DataOutputStream dos, Data d) throws IOException {
			dos.writeInt(d.pair_id);
			dos.writeInt(d.marker_id);
			for(int i=0;i<d.ibd.length;i++) {
				dos.writeFloat(d.ibd[i]);
				}
			}
		
		
		@Override
		public AbstractDataCodec<Data> clone() {
			return new DataCodec();
			}
		
		}
	
	/** marker */
	private class Marker implements Locatable
		{
		final int tid;
		final int pos;
		final String name;
		int index = -1;

		Marker(int tid,int pos,String name) {
			this.tid = tid;
			this.pos = pos;
			this.name = name;
			}
		public String getName() {
			return name;
			}
		@Override
		public String getContig() {
			return dictionary.getSequence(this.tid).getContig();
			}
		public int getPos() {
			return pos;
			}
		@Override
		public int getStart() {
			return getPos();
			}
		@Override
		public int getEnd() {
			return getPos();
			}
		@Override
		public String toString() {
			return getContig()+":"+getPos()+":"+ getName();
			}
		}
	
	private Optional<Sample> findSampleByFamName(final String family,final String name) {
		return this.all_samples.
				stream().
				filter(S->S.family.equals(family) && S.name.equals(name)).
				findFirst();
		}
	
	private Pair findPairBySample(final Sample sn1,final Sample sn2) {
		final Pair pair = new Pair(sn1, sn2);
		final Pair pair2 = this.pair2pair.get(pair);
		if(pair2!=null) return pair2;
		pair.index = this.all_pairs.size();
		this.all_pairs.add(pair);
		this.pair2pair.put(pair,pair);
		return pair;
		}

	
	private void loadMarkers() throws IOException {
		try(BufferedReader br= IOUtils.openPathForBufferedReading(this.markerBedPath)) {
			String line;
			while((line=br.readLine())!=null) {
				if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
				final String[] tokens = CharSplitter.TAB.split(line);
				if(tokens.length!=4) throw new JvarkitException.TokenErrors(4, tokens);
				final String contig = tokens[0];
				final SAMSequenceRecord ssr = this.dictionary.getSequence(tokens[0]);
				if(ssr==null) throw new JvarkitException.ContigNotFoundInDictionary( contig, this.dictionary);

				int start = Integer.parseInt(tokens[1]);
				if(start<1 || start> ssr.getLengthOnReference()) {
					throw new IndexOutOfBoundsException("pos 0<"+start+"<"+ssr.getSequenceLength()+" in "+line);
					}
				final String name = tokens[3];
				final Marker marker = new Marker(ssr.getSequenceIndex(), start+1, name);
				this.all_markers.add(marker);
				}
			}
		if(this.all_markers.isEmpty()) throw new IOException("Not enough markers in "+this.markerBedPath);
		Collections.sort(this.all_markers,new ContigDictComparator(this.dictionary).createLocatableComparator());
		for(int i=0;i< this.all_markers.size();i++) {
			this.all_markers.get(i).index = i ;
			}
		LOG.info("markers N="+this.all_markers.size());
		}
	
	private void loadPedigree() throws IOException {
		try(BufferedReader br= IOUtils.openPathForBufferedReading(this.famPath)) {
			String line;
			while((line=br.readLine())!=null) {
				if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
				final String[] tokens = SPACES.split(line);
				if(tokens.length<2) throw new JvarkitException.TokenErrors(2, tokens);
				if(findSampleByFamName(tokens[0], tokens[1]).isPresent()) {
					throw new IOException("sample was found twice "+line);
					}
				final Sample sample = new Sample(tokens[0], tokens[1]);
				this.all_samples.add(sample);
				}
			}
		if(this.all_samples.size()<2) throw new IOException("Not enough samples in "+this.famPath);
		LOG.info("samples N="+this.all_samples.size());
		}
	
	private void loadIbd(final SortingCollection<Data> sorting,final Path ibdFilePath) throws IOException {
		try(BufferedReader br= IOUtils.openPathForBufferedReading(ibdFilePath)) {
			String line = br.readLine();
			if(line==null) throw new IOException("First line missing in "+ ibdFilePath );
			final String tokens[] = CharSplitter.TAB.split(line);
			if(tokens.length<5) throw new JvarkitException.TokenErrors(5, tokens);
			int num_pairs = Integer.parseInt(tokens[1]);
			final String contig = tokens[2];
			final int num_markers = Integer.parseInt(tokens[3]);
			final Map<String, Marker> name2marker = new HashMap<>(num_markers);
			this.all_markers.stream().
						filter(M->M.getContig().equals(contig)).
						forEach(M->name2marker.put(M.getName(), M));
			
			final List<Marker> markers = new ArrayList<>(num_markers);
			for(int i=4; i< tokens.length;i++) {
				final Marker marker = name2marker.get(tokens[i]);
				if(marker==null) throw new IOException("Marker "+tokens[i]+" not found on chromosome "+contig+" in "+this.markerBedPath);
				markers.add(marker);
				}
			if(markers.size()!=num_markers) {
				throw new IOException("Expected "+num_markers+" in "+ ibdFilePath+" but got "+markers.size());
				}
			
			while((line=br.readLine())!=null) {
				num_pairs--;
				final String[] tokens2 = SPACES.split(line);
				if(tokens2.length!=4+num_markers*3) throw new IOException("Expected at least "+(4+num_markers*3)+" tokens in "+line+" but got "+tokens.length);
				final Sample sample1 = findSampleByFamName(tokens2[0],tokens2[1]).orElseThrow(()->new IllegalArgumentException("sample not found "+tokens2[0]+"/"+tokens2[1]));;
				final Sample sample2 = findSampleByFamName(tokens2[2],tokens2[3]).orElseThrow(()->new IllegalArgumentException("sample not found "+tokens2[2]+"/"+tokens2[3]));
				final Pair pair = findPairBySample(sample1, sample2);
				
				for(int j=0;j< num_markers;j++) {
					final Data data = new Data();
					data.pair_id = pair.index;
					data.marker_id = markers.get(j).index;
					for(int k=0;k<3;++k) {
						data.ibd[k] = Float.parseFloat(tokens2[4+j*3+k]);
					}
					sorting.add(data);
				}
			}
			if(num_pairs!=0) throw new IOException("bad number pairs vs expected");
		}
	}
	
	
	@Override
	public int doWork(final List<String> args) {
		SortingCollection<Data> sorting = null;
		try {
			final List<Path> inputs = IOUtils.unrollPaths(args);
			if(inputs.isEmpty()) {
				LOG.info("IBD file(s) missing");
				return -1;
				}
			this.dictionary =SequenceDictionaryUtils.extractRequired(this.dictPath);
			this.loadMarkers();
			this.loadPedigree();
			
			
			
			sorting = SortingCollection.newInstance(
					Data.class,new DataCodec(),
					(A,B)->A.compare1(B),
					writingSortingCollection.getMaxRecordsInRam(),
					writingSortingCollection.getTmpPaths()
					);
			for(final Path ibdPath : inputs) {
				loadIbd(sorting, ibdPath);
			}

			
			sorting.doneAdding();
			sorting.setDestructiveIteration(true);
			
			
			final Set<VCFHeaderLine> metaData = new HashSet<>();
			VCFFormatHeaderLine ibdHeaderLine = new VCFFormatHeaderLine("I", 1, VCFHeaderLineType.Float,"IBD Value 0");
			metaData.add(ibdHeaderLine);
			ibdHeaderLine = new VCFFormatHeaderLine("J", 1, VCFHeaderLineType.Float,"IBD Value 1");
			metaData.add(ibdHeaderLine);
			ibdHeaderLine = new VCFFormatHeaderLine("K", 1, VCFHeaderLineType.Float,"IBD Value 2");
			metaData.add(ibdHeaderLine);
			final VCFHeader header= new VCFHeader(metaData, 
					this.all_pairs.stream().map(P->P.toString()).collect(Collectors.toList()));
			header.setSequenceDictionary(this.dictionary);
			JVarkitVersion.getInstance().addMetaData(this, header);
			try(VariantContextWriter vcw = super.openVariantContextWriter(null)) {
				vcw.writeHeader(header);
				final List<Allele> alleles  = Collections.singletonList(Allele.REF_N);
				try(CloseableIterator<Data> iter  = sorting.iterator() ) {
					final EqualRangeIterator<Data> eqiter = new EqualRangeIterator<>(iter, (A,B)->A.compare2(B));
					while(eqiter.hasNext()) {
						final List<Data> row = eqiter.next();
						final Marker marker = this.all_markers.get(row.get(0).marker_id);
						final List<Genotype> genotypes = new ArrayList<>(row.size());
						for(int i=0;i< row.size();i++) {
							final Pair pair = this.all_pairs.get(row.get(i).pair_id);
							final GenotypeBuilder gb= new GenotypeBuilder(pair.getName());
							gb.attribute("I", row.get(i).ibd[0]);
							gb.attribute("J", row.get(i).ibd[1]);
							gb.attribute("K", row.get(i).ibd[2]);
							genotypes.add(gb.make());
							}
						final VariantContextBuilder vcb = new VariantContextBuilder(
								null,marker.getContig(),marker.getPos(),marker.getPos(),
								alleles
							).genotypes(genotypes).id(marker.getName());
						vcw.add(vcb.make());
						}
					eqiter.close();
					}
				}
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		
		}
	
	public static void main(final String[] args) {
		new IbdToVcf().instanceMainWithExit(args);
		}
	}
