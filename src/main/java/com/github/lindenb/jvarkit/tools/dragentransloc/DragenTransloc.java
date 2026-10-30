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
package com.github.lindenb.jvarkit.tools.dragentransloc;

import java.io.IOException;
import java.nio.file.Path;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.SmartComparator;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.variant.variantcontext.Breakend;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC

# Motivation 

find translocations using Dragen , filter out those found in another set of VCFS (controls)


e.g. true translocation found in dragen:

```
chrA  56120849  DRAGEN:BND:33583:0:1:0:12:0:1  G  G[chrB:174700[    999  PASS  SVTYPE=BND;MATEID=DRAGEN:BND:33583:0:1:0:12:0:0;
chrA  56120857  DRAGEN:BND:33583:0:1:1:7:0:1   G  ]chrB:174678]G    999  PASS  SVTYPE=BND;MATEID=DRAGEN:BND:33583:0:1:1:7:0:0;
chrB  174678    DRAGEN:BND:33583:0:1:1:7:0:0   T  T[chrA:56120857[  999  PASS  SVTYPE=BND;MATEID=DRAGEN:BND:33583:0:1:1:7:0:1;
chrB  174700    DRAGEN:BND:33583:0:1:0:12:0:0  C  ]chrA:56120849]C  999  PASS  SVTYPE=BND;MATEID=DRAGEN:BND:33583:0:1:0:12:0:1;
```

# Example

```
mkdir TMP

find /path -type f -name "*.sv.vcf.gz" > controls.list

java -jar jvarkit.jar dragentransloc \\
	--controls controls.list \\
	--distance 100 \\
	"in.vcf" > out.vcf

```

END_DOC
 */
@Program(
	name="dragentransloc",
	description="Finds dragen translocations",
	keywords={"vcf","dragen","inversion","bnd","sv"},
	creationDate="20261028",
	modificationDate="20261028")
public class DragenTransloc extends OnePassVcfLauncher {
	private static final String MATEID ="MATEID";
	private static Logger LOG=Logger.of(DragenTransloc.class);
	@Parameter(names={"--control","--controls"},description="another vcf or a text file with the '.list' suffix containing the path to the controls. Any translocation will be removed if it's found in a control ")
	private Path control_source =null;
	@Parameter(names={"--distance"},description="two BND on the same chromosome at a distance of 'x' are considered the same. " + DistanceParser.OPT_DESCRIPTION,converter = DistanceParser.StringConverter.class,splitter = NoSplitter.class)
	private int bnd_distance=10;
	@Parameter(names={"--blacklist-bed","--exclude-bed"},description="Do not use the variant falling in that region. ")
	private Path exclude_bed_path = null;
	@Parameter(names={"--accept-filtered"},description="Accept FILTERed variant")
	private boolean accept_filtered_variants=false;
	@Parameter(names={"--debug"},description="debug")
	private boolean debug_flag=false;
	@Parameter(names={"--disable-discordant-contig","--disable-same-contig","--disable-same-chromosome"},description="accept BND where both break-end are on the same chromosome.")
	private boolean accept_same_contigs_flags=false;
	@Parameter(names={"--min-dp","--min-depth"},description="min value for INFO/BND_DEPTH and INFO/MATE_BND_DEPTH. -1 = ignore")
	private int min_bnd_depth=-1;
	@Parameter(names={"--skip-imprecise"},description="Skip variants with the INFO/IMPRECISE flag.")
	private boolean skip_imprecise_flag = false;
	@Parameter(names={"--two-junctions-factor"},description="When looking for the pair (ID/MATEID) in the main vcf, there should be another pair of  (ID/MATEID) for the other side of the breakend, mapping on the same chromosome and at a (--distance)*'x'."
			+ "Bute one of the two junctions may be removed when found in a control VCF. Ignored if <1.")
	private int two_junctions_factor = -1;

	
	private final List<Path> controls_path=new ArrayList<>();
	
	private static class Translocation extends AbstractList<VariantContext> {
		final List<VariantContext> variants= new ArrayList<>(2);
		@Override
		public VariantContext get(int index) {
			return variants.get(index);
			}
		@Override
		public int size() {
			return variants.size();
			}
		
		boolean hasACommonBreakEndWith(final Translocation other,final int distance) {
			for(VariantContext vc1: this.variants) {
				final Locatable extended1 = extendToCIPOS(vc1);
				for(VariantContext vc2: other.variants) {
					final Locatable extended2 = extendToCIPOS(vc2);
					if(extended1.contigsMatch(extended2) && extended1.withinDistanceOf(extended2, distance)) {
						return true;
						}
					}
				}
			return false;
			}

		boolean sameAs(final Translocation other,final int distance) {
			for(VariantContext vc1: this.variants) {
				boolean ok=false;
				final Locatable extended1 = extendToCIPOS(vc1);
				for(VariantContext vc2: other.variants) {
					final Locatable extended2 = extendToCIPOS(vc2);
					if(extended1.contigsMatch(extended2) && extended1.withinDistanceOf(extended2, distance)) {
						ok=true;
						break;
						}
					}
				if(!ok) return false;
				}
			return true;
			}
		@Override
		public String toString() {
			return "Translocation "+variants.stream().map(V->V.getContig()+":"+V.getStart()).collect(Collectors.joining(" / "));
			}
		}
	/* use CIPOS to extend interval */
	private static Locatable extendToCIPOS(final VariantContext ctx) {
		if(!ctx.hasAttribute("CIPOS")) return ctx;
		final int extend = ctx.getAttributeAsIntList("CIPOS", 0).stream().mapToInt(I->Math.abs(I.intValue())).max().orElse(0);
		
		return new SimpleInterval(ctx.getContig(),
				Math.max(1, ctx.getStart()-extend),
				ctx.getEnd()+extend
				);
		}
	
	
	/** extract translocation from a vcf file */
	private List<Translocation> extract(final CloseableIterator<VariantContext> iter,
			boolean accept_filtered,
			final IntervalTreeMap<Boolean> limit
			) {
		final Map<String,Translocation> id2transloc= new HashMap<>();
		while(iter.hasNext()) {
			final VariantContext ctx = iter.next();
			
			
			if(!accept_filtered && ctx.isFiltered()) {
				continue;
				}
			if(!ctx.hasID()) continue;
			if(!ctx.hasAttribute(VCFConstants.SVTYPE)){
				continue;
				}
			if(!ctx.hasAttribute(MATEID)) {
				continue;
				}
			
			if(!ctx.getAttributeAsString(VCFConstants.SVTYPE, "").equals("BND")) {
				continue;
				}
			
			/*
			if(ctx.hasGenotypes()) {
				final Genotype gt = ctx.getGenotype(0);
				if(gt.hasAnyAttribute("VF") && gt.isHet()) {
					final int[] L=  CharSplitter.COMMA.stream(String.valueOf(gt.getAttributeAsString("VF","")))
							.mapToInt(S->Integer.valueOf(S)).toArray();
					if(L.length==2) {
						double f= (L[0]*1.0)/(L[0]+L[1]);
						final double t=0.25;
						if(f < t || f>(1.0-t)) continue;
						}
					}
				}*/
			
			if(this.skip_imprecise_flag && ctx.hasAttribute("IMPRECISE")) {
				continue;
				}
			
			if(this.min_bnd_depth > -1) {
				if(ctx.hasAttribute("BND_DEPTH")) {
					final int dp = ctx.getAttributeAsInt("BND_DEPTH", 0);
					if(dp< this.min_bnd_depth) continue;
					}
				if(ctx.hasAttribute("MATE_BND_DEPTH")) {
					final int dp = ctx.getAttributeAsInt("MATE_BND_DEPTH", 0);
					if(dp< this.min_bnd_depth) continue;
					}
				}
			
			if(limit!=null && !limit.containsOverlapping(ctx)) {
				continue;
				}
			
			final Optional<Breakend> optBnd =  Breakend.parse(ctx.getAlleles().get(1));
			if(!optBnd.isPresent()){
				continue;
				}
			final Breakend bnd = optBnd.get();
			if(!this.accept_same_contigs_flags && bnd.contigsMatch(ctx)) {
				continue;
				}
			final String mateid = ctx.getAttributeAsString(MATEID, "");
			if(StringUtils.isBlank(mateid)) {
				continue;
				}
			Translocation tr = id2transloc.get(mateid);
			if(tr==null) {
				tr=new Translocation();
				id2transloc.put(ctx.getID(), tr);
				}
			tr.variants.add(ctx);
			}
		final List<Translocation> L= new ArrayList<>(id2transloc.values());
		L.removeIf(TR->TR.size()!=2);
		return L;
		}
	
	
	@Override
	protected int beforeVcf() {
		if(this.control_source!=null) {
			this.controls_path.addAll( IOUtils.unrollPath(this.control_source));
			}
		return super.beforeVcf();
		}
	
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}	
	
	private void applyTwoJunctionsRequired(final List<Translocation> translocs) {
		final  List<Translocation> L=new ArrayList<>(translocs.size());
		for(int x=0;x< translocs.size();++x) {
			for(int y=0;y< translocs.size();++y) {
				if(x==y) continue;
				if( translocs.get(x).sameAs(translocs.get(y), this.bnd_distance*this.two_junctions_factor)) {
					L.add(translocs.get(x));
					break;
					}
				}
			}
		translocs.clear();
		translocs.addAll(L);
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iter, VariantContextWriter out) {
		final VCFHeader header= iter.getHeader();
		final SAMSequenceDictionary dict1 = header.getSequenceDictionary(); // may be null
		
		
		if(header.getInfoHeaderLine(MATEID)==null) {
			LOG.error("INFO/"+MATEID+" is not defined in "+inputName);
			return -1;
			}
		
		if(header.getNGenotypeSamples()>1) {
			LOG.error("Cannot handle more than one genotyped sample");
			return -1;
			}
		final List<Translocation> translocs = extract(iter,this.accept_filtered_variants,null);
		if(translocs.isEmpty()) {
			LOG.warn("no translocation was found in input");
			}
		
		
		/* apply exclude map */
		if(this.exclude_bed_path!=null && !translocs.isEmpty()) {
			try(BedLineReader blr= BedLineReader.openPath(this.exclude_bed_path)) {
				if(dict1!=null) blr.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict1));
				final IntervalTreeMap<Boolean> excludeMap = blr.toBooleanIntervalTreeMap();
				translocs.removeIf(TR->TR.variants.stream().anyMatch(VC->excludeMap.containsOverlapping(VC)));
				}
			}
		
		/* check we have a junction for each side of the translocation */
		if(this.two_junctions_factor>=1) {
			applyTwoJunctionsRequired(translocs);
			}

		
		
		if(this.debug_flag) {
			for(Translocation tr: translocs) {
				LOG.debug("found case transloc "+tr);
				}
			}
		
		IntervalTreeMap<Boolean> limitOverlap = null;
		
		boolean multiple_dict_flag=false;
		for(Path ctrl_path : this.controls_path) {
			if(translocs.isEmpty()) break;
			
			try(VCFReader vcfReader= new VCFFileReader(ctrl_path,false)) {
				final VCFHeader header2= vcfReader.getHeader();
				if(header2.getInfoHeaderLine(MATEID)==null) {
					LOG.error("INFO/"+MATEID+" is not defined in "+ctrl_path+". skipping.");
					continue;
					}
				final SAMSequenceDictionary dict2 = header2.getSequenceDictionary();
				if(dict1!=null && dict2!=null) {
					if(!SequenceUtil.areSequenceDictionariesEqual(dict1, dict2) && !multiple_dict_flag) {
						multiple_dict_flag = true;
						LOG.warn("NOT the same dictionaries between "+inputName +" and "+ctrl_path+" !. Use with caution. This message is printed one time.");
						}
					}
				
				// create an intervaltreemap to limit the number of tests
				if(limitOverlap==null) {
					limitOverlap = new IntervalTreeMap<>();
					for(Translocation tr:translocs) {
						for(VariantContext ctx: tr) {
							final Locatable extended = extendToCIPOS(ctx);
							limitOverlap.put(
								new Interval(
										extended.getContig(),
										Math.max(1, extended.getStart()- 2*this.bnd_distance),
										extended.getEnd()+2*this.bnd_distance
										)	
								, Boolean.TRUE);
							}
						}
					}
				
				final List<Translocation> control_translocs;
				try(CloseableIterator<VariantContext> iter2 = vcfReader.iterator()) {
					control_translocs =  extract(iter2,true /* always accept FILTERed variant in controls*/,limitOverlap);
					}
				if(control_translocs.isEmpty()) continue;
				// scan all transloc in the case, remove if any has a match in the control
				int x=0;
				while(x< translocs.size()) {
					final Translocation tr= translocs.get(x);
					final Translocation same= control_translocs.stream()
							.filter(CT->CT.hasACommonBreakEndWith(tr,bnd_distance))
							.findAny()
							.orElse(null);
					
					if(same!=null) {
						if(this.debug_flag) {
							LOG.debug("removing "+tr+" because it's the same than "+same+" in "+ctrl_path+" remains "+(translocs.size()-1));
							}
						//force re-build overlap next time
						limitOverlap = null;
						translocs.remove(x);
						}
					else
						{
						x++;
						}
					}
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		
		/* check we have a junction for each side of the translocation */
		if(this.two_junctions_factor>=1) {
			applyTwoJunctionsRequired(translocs);
			}
		
		/* create VariantContextComparator to sort the output variants */
		final VariantContextComparator cmp;
		if(dict1!=null) {
			cmp = new VariantContextComparator(dict1);
			}
		else if(translocs.isEmpty()) {
			cmp = new VariantContextComparator(Collections.singletonList("dummy_contig"));
			}
		else
			{
			cmp = new VariantContextComparator(
					translocs.stream()
						.flatMap(TR->TR.variants.stream())
						.map(V->V.getContig())
						.collect(Collectors.toSet())
						.stream()
						.sorted(new SmartComparator())
						.collect(Collectors.toList())
					);
			}
		/* write variants */
		JVarkitVersion.getInstance().addMetaData(this, header);
		out.writeHeader(header);
		for(VariantContext ctx: translocs.stream()
				.flatMap(TR->TR.variants.stream())
				.sorted(cmp)
				.collect(Collectors.toList())) {
			out.add(ctx);
			}
		return 0;
		}
	
	
public static void main(final String[] args) {
	new DragenTransloc().instanceMainWithExit(args);
	}
}
