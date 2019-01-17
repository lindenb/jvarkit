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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;

@Program(name="vcfburden",
	description="Solena: vcf to (chrom/pos/ref/alt/individus(G:0/1/2/-9)",
	deprecatedMsg="deprecated")
public class VcfBurden extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfBurden.class).make();


	@Parameter(names={"-o","--output"},description="zip file")
	private File outputFile = null;
	@Parameter(names="-H",description=" only high damage")
	private boolean highdamage=false;
	private Map<String,Boolean> _gene2seen=null;
	@Parameter(names="-d",description=" (dir) base zip dir")
	private String dirName="burden";
	@Parameter(names="-f",description=" print ALL consequence+ colum, SO terms (mail matilde 11/10/2015 11:38 AM)")
	private boolean printSOTerms = false;//mail matilde 11/10/2015 11:38 AM
	@Parameter(names="-p",description=" print position in CDS")
	private boolean printPositionInCDS = false;
	@Parameter(names="-q",description=" print VQSLOD")
	private boolean printVQSLOD=false;
	
	
	@Parameter(names="-g",description=" (file) optional list of gene names (restrict genes, print genes without data)")
	private File userlistOfGene= null;
	
	
	

	
	
	public VcfBurden()
		{
		}
	
	private static class VariantAndCsq
		{
		final VariantContext variant;
		final Set<SequenceOntologyTree.Term> terms;
		final Integer pos_in_cds;
		final Float vqslod;
		VariantAndCsq(VariantContext variant,final Set<SequenceOntologyTree.Term> terms,final Integer pos_in_cds, final Float vqslod)
			{
			this.variant = variant;
			this.terms = new HashSet<>(terms);
			this.pos_in_cds=pos_in_cds;
			this.vqslod = vqslod;
			}
		
		Integer getPositionInCDS()
			{
			return this.pos_in_cds;
			}
		@Override
		public String toString() {
			return variant.getContig()+":"+variant.getStart()+"-"+variant.getEnd()+" "+
					 variant.getReference()+"/"+variant.getAlternateAlleles();
			}
		}
	
	private final Comparator<VariantAndCsq> variantAndCsqComparator = new Comparator<VariantAndCsq>()
		{
		@Override
		public int compare(
				final VariantAndCsq v1,
				final VariantAndCsq v2) {
			final VariantContext o1 = v1.variant;
			final VariantContext o2 = v2.variant;
			int i = o1.getContig().compareTo(o2.getContig());
			if(i!=0) return i;
			i = o1.getStart() - o2.getStart();
			if(i!=0) return i;
			i = o1.getReference().compareTo(o2.getReference());
			if(i!=0) return i;
			i = o1.getAlternateAllele(0).compareTo(o2.getAlternateAllele(0));
			if(i!=0) return i;
			return 0;
			}
		};
	
	
	private static class GeneTranscript
		{
		String geneName;
		String transcriptName;
		GeneTranscript(String geneName,String transcriptName)
			{
			this.geneName = geneName;
			this.transcriptName = transcriptName;
			}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result
					+ ((geneName == null) ? 0 : geneName.hashCode());
			result = prime * result
					+ ((transcriptName == null) ? 0 : transcriptName.hashCode());
			return result;
			}
		@Override
		public boolean equals(Object obj) {
			if (this == obj) return true;
			if (obj == null) return false;
			if (getClass() != obj.getClass()) return false;
			GeneTranscript other = (GeneTranscript) obj;
			if (geneName == null) {
				if (other.geneName != null)
					return false;
			} else if (!geneName.equals(other.geneName))
				return false;
			if (transcriptName == null) {
				if (other.transcriptName != null)
					return false;
			} else if (!transcriptName.equals(other.transcriptName))
				return false;
			return true;
			}
		
		}

	
	
	private void dumpVariants(
			final ZipOutputStream zout,
			final String contig,
			final String filename,
			final List<String> samples,
			final List<VariantAndCsq> variants
			) throws IOException
		{
		for(int i=0;i<variants.size();++i ) {
			for(int j=i+1;j<variants.size();++j ) {
				VariantAndCsq v1  = variants.get(i);
				VariantAndCsq v2  = variants.get(j);
				if(this.variantAndCsqComparator.compare(v1, v2)==0) {
					throw new IOException(v1+" "+v2);
				}
			}
		}
		
		
		int outCount=0;
		LOG.info(this.dirName+"/"+contig+"_"+filename+".txt");
		final ZipEntry ze = new ZipEntry(this.dirName+"/"+contig+"_"+filename+".txt");
		zout.putNextEntry(ze);
		final  PrintWriter pw = new PrintWriter(zout);
		pw.print("CHROM\tPOS\tREF\tALT");
		
		if(printVQSLOD)
			{
			pw.print("\t");
			pw.print("VQSLOD");
			}
		
		if(printPositionInCDS)
			{
			pw.print("\t");
			pw.print("PositionInCDS");
			}
		
		if(printSOTerms)
			{
			pw.print("\t");
			pw.print("SO");//TODO
			}
		for(final String sample:samples)
			{
			pw.print("\t");
			pw.print(sample);
			}
		pw.println();
		for(final VariantAndCsq vac:variants)
			{
			pw.print(vac.variant.getContig());
			pw.print("\t");
			pw.print(vac.variant.getStart());
			pw.print("\t");
			pw.print(vac.variant.getReference().getDisplayString());
			pw.print("\t");
			pw.print(vac.variant.getAlternateAlleles().get(0).getDisplayString());
			
			
			if(printVQSLOD)
				{
				pw.print("\t");
				pw.print(String.valueOf(vac.vqslod));
				}

			
			if(printPositionInCDS)
				{
				pw.print("\t");
				pw.print(vac.getPositionInCDS()==null?"N/A":""+vac.getPositionInCDS());
				}
			
			if(printSOTerms)
				{
				pw.print("\t");
				for(SequenceOntologyTree.Term t:vac.terms)
					{
					pw.print(t.getAcn());
					pw.print(";");
					}
				}
			for(final String sample:samples)
				{
				final Genotype g=vac.variant.getGenotype(sample);
				pw.print("\t");
				if(g.isHomRef())
					{
					pw.print("0");
					}
				else if(g.isHomVar())
					{
					pw.print("2");
					}
				else if(g.isHet())
					{
					pw.print("1");
					}
				else
					{
					pw.print("-9");
					}
				}
			pw.println();
			outCount++;
			}
		pw.flush();
		zout.closeEntry();
		LOG.info(this.dirName+"/"+contig+"_"+filename+".txt N="+outCount);
		}
	
	private int doWork2(List<String> args) {
		ZipOutputStream zout=null;
		FileOutputStream fout=null;
		VCFIterator in=null;
		try
			{
			if(args.isEmpty())
				{
				LOG.info("reading from stdin.");
				in = VCFUtils.createVCFIteratorStdin();
				}
			else if(args.size()==1)
				{
				String filename=args.get(0);
				LOG.info("reading from "+filename);
				in = VCFUtils.createVCFIterator(filename);
				}
			else
				{
				LOG.error("Illegal number of arguments.");
				return -1;
				}
			if(outputFile==null)
				{
				LOG.error("undefined output");
				return -1;
				}
			else if(!outputFile.getName().endsWith(".zip"))
				{
				LOG.error("output "+outputFile+" should end with .zip");
				return -1;
				}
			else
				{
				fout = new FileOutputStream(outputFile);
				zout = new ZipOutputStream(fout);
				}
			final List<String> samples= in.getHeader().getSampleNamesInOrder();
			final VCFHeader header=in.getHeader();
			String prev_chrom = null;
			final VepPredictionParser vepPredParser=new VepPredictionParserFactory().header(header).get();
			final Map<GeneTranscript,List<VariantAndCsq>> gene2variants=new HashMap<>();
			final SequenceOntologyTree soTree= SequenceOntologyTree.getInstance();
			final Set<SequenceOntologyTree.Term> acn=new HashSet<>();
			/* mail solena  *SO en remplacement des SO actuels (VEP HIGH + MODERATE) - pas la peine de faire retourner les analyses mais servira pour les futures analyses burden */
			String acn_list[]=new String[]{
					"SO:0001893",  "SO:0001574",  "SO:0001575", 
					"SO:0001587",  "SO:0001589",  "SO:0001578", 
					"SO:0002012",  "SO:0001889",  "SO:0001821", 
					"SO:0001822",  "SO:0001583",  "SO:0001818"
					
					/*
					"SO:0001589", "SO:0001587", "SO:0001582", "SO:0001583",
					"SO:0001575", "SO:0001578", "SO:0001574", "SO:0001889",
					"SO:0001821", "SO:0001822", "SO:0001893"*/
					};
			
			if(this.highdamage)
				{
				acn_list=new String[]{
						"SO:0001893",  "SO:0001574",  "SO:0001575",
						"SO:0001587",  "SO:0001589",  "SO:0001578", 
						"SO:0002012",  "SO:0001889"
						};
				}
			
			for(final String acns:acn_list)
				{
				final SequenceOntologyTree.Term  tacn = soTree.getTermByAcn(acns);
				if(tacn==null)
					{
					in.close();
					throw new NullPointerException("tacn == null pour "+acns);
					}
				acn.addAll(tacn.getAllDescendants());
				}
			
		
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(in.getHeader());
			for(;;)
				{
				VariantContext ctx1=null;
				if(in.hasNext())
					{
					ctx1= progress.watch(in.next());
					if(ctx1.getAlternateAlleles().size()!=1)
						{
						//info("count(ALT)!=1 in "+ctx1.getChr()+":"+ctx1.getStart());
						continue;
						}
					if(ctx1.getAlternateAlleles().get(0).equals(Allele.SPAN_DEL)) {
						continue;
						}
					}
				
				if(ctx1==null || !ctx1.getContig().equals(prev_chrom))
					{
					LOG.info("DUMP to zip n="+gene2variants.size());
					final Set<String> geneNames= new HashSet<>();
					for(GeneTranscript gene_transcript:gene2variants.keySet() )
						{
						geneNames.add(gene_transcript.geneName);
						
						final Set<VariantAndCsq> uniq =new TreeSet<>(this.variantAndCsqComparator);
						uniq.addAll(gene2variants.get(gene_transcript));
						
						dumpVariants(zout,
								prev_chrom,
								gene_transcript.geneName+"_"+gene_transcript.transcriptName,
								samples,
								new ArrayList<VariantAndCsq>(uniq)
								);
						LOG.info("dumped" +gene_transcript.geneName);
						}
					LOG.info("loop over geneName");
					for(final String geneName : geneNames)
						{
						
						final SortedSet<VariantAndCsq> lis_nm = new TreeSet<>(this.variantAndCsqComparator);
						final SortedSet<VariantAndCsq> lis_all = new TreeSet<>(this.variantAndCsqComparator);
						final SortedSet<VariantAndCsq> lis_refseq = new TreeSet<>(this.variantAndCsqComparator);
						final SortedSet<VariantAndCsq> lis_enst = new TreeSet<>(this.variantAndCsqComparator);
						LOG.info("loop over gene2variants");
						for(final GeneTranscript gene_transcript:gene2variants.keySet() )
							{
							if(!geneName.equals(gene_transcript.geneName)) continue;
							lis_all.addAll(gene2variants.get(gene_transcript));
							if(gene_transcript.transcriptName.startsWith("NM_"))
								{
								lis_nm.addAll(gene2variants.get(gene_transcript));
								}
							if(! gene_transcript.transcriptName.startsWith("ENST"))
								{
								lis_refseq.addAll(gene2variants.get(gene_transcript));
								}
							if( gene_transcript.transcriptName.startsWith("ENST"))
								{
								lis_enst.addAll(gene2variants.get(gene_transcript));
								}
							}
						LOG.info("dump_ALL_TRANSCRIPTS");
						dumpVariants(zout,
								prev_chrom,
								geneName+"_ALL_TRANSCRIPTS",
								samples,
								new ArrayList<VariantAndCsq>(lis_all)
								);
						LOG.info("dump_ALL_NM");
						dumpVariants(zout,
								prev_chrom,
								geneName+"_ALL_NM",
								samples,
								new ArrayList<VariantAndCsq>(lis_nm)
								);
						LOG.info("dump_ALL_REFSEQ");
						dumpVariants(zout,
								prev_chrom,
								geneName+"_ALL_REFSEQ",
								samples,
								new ArrayList<VariantAndCsq>(lis_refseq)
								);
						LOG.info("dump_ALL_ENST");
						dumpVariants(zout,
								prev_chrom,
								geneName+"_ALL_ENST",
								samples,
								new ArrayList<VariantAndCsq>(lis_enst)
								);
						}
					
					if(ctx1==null) break;
					LOG.info("gene2variants");
					gene2variants.clear();
					LOG.info("System.gc();");
					System.gc();
					prev_chrom = ctx1.getContig();
					LOG.info("prev_chrom="+prev_chrom);
					}
				final Set<GeneTranscript> seen_names=new HashSet<>();
				for(final VepPredictionParser.VepPrediction pred: vepPredParser.getPredictions(ctx1))
					{
					String geneName= pred.getSymbol();
					if(geneName==null || geneName.trim().isEmpty()) continue;

					
					if(this._gene2seen!=null)
						{
						if(!this._gene2seen.containsKey(geneName)) continue;
						
						}
					final String transcriptName = pred.getFeature();
					if(transcriptName==null || transcriptName.trim().isEmpty())
						{
						LOG.info("No transcript in "+ctx1);
						continue;
						}
					
					final GeneTranscript geneTranscript = new GeneTranscript(geneName, transcriptName);
					
					if(seen_names.contains(geneTranscript)) continue;
					boolean ok=false;
					for(SequenceOntologyTree.Term so:pred.getSOTerms())
						{
						if(acn.contains(so))
							{
							ok=true;
							}
						}
					if(!printSOTerms && !ok) continue;
					
					List<VariantAndCsq> L = gene2variants.get(geneTranscript);
					if(L==null)
						{
						L=new ArrayList<>();
						gene2variants.put(geneTranscript,L);
						}
					Float vqslod = null;
					if( this.printVQSLOD && ctx1.hasAttribute("VQSLOD")){
						vqslod = (float)ctx1.getAttributeAsDouble("VQSLOD",-9999999.0);	
					}
					
					L.add(new VariantAndCsq(
							ctx1,
							pred.getSOTerms(),
							pred.getPositionInCDS(),
							vqslod
							));
					seen_names.add(geneTranscript);
					if(this._gene2seen!=null)
						{
						this._gene2seen.put(geneTranscript.geneName, Boolean.TRUE);
						}
					}
				
				}
			
			if(this._gene2seen!=null)
				{
				final List<VariantAndCsq> emptylist = Collections.emptyList();
				for(final String gene:this._gene2seen.keySet())
					{
					if(this._gene2seen.get(gene).equals(Boolean.TRUE)) continue;
					LOG.warning("Gene not found : "+gene);
					dumpVariants(zout,
							"UNDEFINED",
							gene+"_000000000000000.txt",
							samples,
							emptylist
							);
					}
				}
			
			progress.finish();
			
			zout.finish();
			fout.flush();
			zout.flush();
			
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			CloserUtil.close(zout);
			CloserUtil.close(fout);
			}
		}
	

	
	
	@Override
	public int doWork(List<String> args)
		{
		
		BufferedReader in=null;
		try {
			if(userlistOfGene!=null) {
				if(this._gene2seen==null) this._gene2seen=new HashMap<>();
				in = IOUtils.openFileForBufferedReading(userlistOfGene);
				String line;
				while((line=in.readLine())!=null)
					{
					line=line.trim();
					if(line.isEmpty()||line.startsWith("#")) continue;
					this._gene2seen.put(line, Boolean.FALSE);
					}
				in.close();
				}
			return doWork2(args);
			}
		catch (Exception e) {
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			}
		}

	public static void main(String[] args)
		{
		new VcfBurden().instanceMainWithExit(args);
		}
	}
