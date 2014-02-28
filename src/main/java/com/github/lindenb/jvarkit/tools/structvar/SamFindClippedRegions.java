package com.github.lindenb.jvarkit.tools.structvar;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broad.tribble.readers.LineIterator;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFFormatHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFHeaderVersion;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.variant.vcf.VCFSimpleHeaderLine;


import net.sf.picard.PicardException;
import net.sf.picard.sam.MergingSamRecordIterator;
import net.sf.picard.sam.SamFileHeaderMerger;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.SequenceUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamSequenceRecordTreeMap;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

public class SamFindClippedRegions extends AbstractCommandLineProgram
	{
	private int min_clip_length=10;
	private boolean ignore_poly_x=false;
	private int rounding_pos=5;
	//private float min_fraction_of_clipped_records=0.3f;
	private int min_depth=10;
	
	private static class SuspectRgn
		{
		int tid=-1;
		int pos=0;
		Map<String,int[]> sample2count=new HashMap<>(); 
		
		SuspectRgn()
			{
			}
		@Override
		public String toString() {
			return "tid="+tid+":"+pos+" N="+sample2count.size();
			}
		}
	
	private class Input
		implements Closeable
		{
		File bamFile;
		SAMFileReader samFileReaderScan;
		SAMFileReader samFileReaderDepth;
		SAMFileHeader header;
		String sampleName=null;
		Input(File bamFile)
			{
			info("Reading from "+bamFile);
			this.bamFile=bamFile;
			this.samFileReaderScan=new SAMFileReader(this.bamFile);
			this.samFileReaderScan.setValidationStringency(ValidationStringency.SILENT);
			this.header=this.samFileReaderScan.getFileHeader();
			for(SAMReadGroupRecord g:this.header.getReadGroups())
				{
				if(g.getSample()==null || g.getSample().isEmpty())
					{
					warning("Read group "+g.getId()+" has no sample in "+bamFile);
					continue;
					}
				if(sampleName!=null && !sampleName.equals(g.getSample()))
					{	
					throw new PicardException("Sorry one bam can only contain one sample and I got "+sampleName+" and "+g.getSample()+" in "+bamFile);
					}
				sampleName=g.getSample();
				}
			if(sampleName==null)
				{
				throw new PicardException("Bam "+this.bamFile+" doesn't contain a sample defined in header/read-group");
				}
			info("Sample is "+sampleName);
			
			this.samFileReaderDepth=new SAMFileReader(this.bamFile);
			this.samFileReaderDepth.setValidationStringency(ValidationStringency.SILENT);	
			}
		
		@Override
		public void close() throws IOException {
			CloserUtil.close(samFileReaderScan);
			CloserUtil.close(samFileReaderDepth);
			}
		}
	
	@Override
	public String getProgramDescription() {
		return "";
		}
	

	
	private boolean isClipped(final SAMRecord rec)
		{
		return isClipped(rec,0) || isClipped(rec,1);
		}
	private CigarElement getCigarElement(final SAMRecord rec,int side)
		{
		Cigar cigar=rec.getCigar();
		if(cigar==null || cigar.isEmpty()) return null;
		CigarElement ce=cigar.getCigarElement(side==0?0:cigar.numCigarElements()-1);
		return ce;
		}
	private boolean isClipped(final SAMRecord rec,int side)
		{
		CigarElement ce=getCigarElement(rec,side);
		if(ce==null) return false;
		return ce.getOperator().equals(CigarOperator.S) && ce.getLength() >= min_clip_length;
		}
	
	private String getClippedSequence(final SAMRecord rec,int side)
		{
		CigarElement ce=getCigarElement(rec,side);
		if(side==0)
			{
			return rec.getReadString().substring(0, ce.getLength());
			}
		else
			{
			return rec.getReadString().substring(rec.getReadLength()-ce.getLength());
			}
		}
	
	
	private static boolean same(char c1,char c2)
		{
		if(c1=='N' || c2=='N') return false;
		return Character.toUpperCase(c1)==Character.toUpperCase(c2);
		}
	
	
	private boolean accept(final SAMRecord rec)
		{
		if(rec.getReadUnmappedFlag()) return false;
		if(rec.isSecondaryOrSupplementary()) return false;
		if(rec.getDuplicateReadFlag()) return false;
		SAMReadGroupRecord g=rec.getReadGroup();
		if(g==null || g.getSample()==null || g.getSample().isEmpty()) return false;
		return true;
		}
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -c (int) min size of clipped read. default:"+min_clip_length);
		out.println(" -x ignore poly-X");
		out.println(" -B (file) bed file (optional)");
		out.println(" -R (ref) "+getMessageBundle("reference.faidx")+". Optional");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		File bedFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"c:xB:"))!=-1)
			{
			switch(c)
				{
				case 'B': bedFile=new File(opt.getOptArg());break;
				case 'x': ignore_poly_x=true;break;
				case 'c':min_clip_length=Integer.parseInt(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
	
		
		if(opt.getOptInd()==args.length)
			{
			error(getMessageBundle("illegal.number.of.arguments"));
			return -1;
			}
		
		List<Input> inputs=new ArrayList<Input>();
		
		VariantContextWriter w=null;
		//SAMFileWriter w=null;
		try
			{
			SAMSequenceDictionary dict=null;

			/* create input, collect sample names */
			Map<String,Input> sample2input=new HashMap<String,Input>();
			for(int optind=opt.getOptInd();
					optind< args.length && inputs.size()<2;//TODO
					++optind)
				{
				Input input=new Input(new File(args[optind]));
				inputs.add(input);
				
				if(sample2input.containsKey(input.sampleName))
					{
					error("Duplicate sample "+input.sampleName+" in "+input.bamFile+" and "+sample2input.get(input.sampleName).bamFile);
					return -1;
					}
				sample2input.put(input.sampleName, input);
				if(dict==null)
					{
					dict=input.header.getSequenceDictionary();
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(dict, input.header.getSequenceDictionary()))
					{
					error("Found more than one dictint sequence dictionary");
					return -1;
					}
					
				}
			info("Sample N= "+sample2input.size());
			
			
			/* create merged iterator */
			List<SAMFileHeader> headers=new ArrayList<SAMFileHeader>(sample2input.size());
			for(Input input:inputs) headers.add(input.header);
			SamFileHeaderMerger headerMerger=new SamFileHeaderMerger(SortOrder.coordinate,headers, true);
			
			List<SAMFileReader> readers=new ArrayList<SAMFileReader>(sample2input.size());
			for(Input input:inputs) readers.add(input.samFileReaderScan);
			
			MergingSamRecordIterator merginIter=new MergingSamRecordIterator(headerMerger, readers, true);
			
			Allele reference_allele=Allele.create("N",true);			
			Allele alternate_allele=Allele.create("<DEL>",false);

			Set<VCFHeaderLine> vcfHeaderLines=new HashSet<VCFHeaderLine>();
			vcfHeaderLines.add(new VCFSimpleHeaderLine(
					"<ID="+alternate_allele.getDisplayString()+",Description=\"StructVar\">",
					VCFHeaderVersion.VCF4_1,
					VCFConstants.ALT_HEADER_START.substring(2),
					Arrays.asList("ID", "Description")
					));
			vcfHeaderLines.add(new VCFInfoHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Approximate read depth."));
			vcfHeaderLines.add(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			vcfHeaderLines.add(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			vcfHeaderLines.add(new VCFFormatHeaderLine("CN",1,VCFHeaderLineType.String,"count(Clipped 5'):count(Clipped 3')"));

			if(dict!=null)
				{
				vcfHeaderLines.addAll(VCFUtils.samSequenceDictToVCFContigHeaderLine(dict));
				}
			
			VCFHeader vcfHeader=new VCFHeader(vcfHeaderLines,sample2input.keySet());
			
			w=VCFUtils.createVariantContextWriterToStdout();
			w.writeHeader(vcfHeader);
			

			
			SamSequenceRecordTreeMap<Boolean> intervals=null;
			//w=swf.make(header, System.out);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
			
			if(bedFile!=null)
				{
				intervals=new SamSequenceRecordTreeMap<Boolean>(dict);
				info("Reading "+bedFile);
				LineIterator r=IOUtils.openFileForLineIterator(bedFile);
				while(r.hasNext())
					{
					String line=r.next();
					if(line.startsWith("#") || line.isEmpty()) continue;
					String tokens[]=line.split("[\t]");
					if(tokens.length<3)
						{
						warning("Bad bed line in "+bedFile+" "+line);
						continue;
						}
					if(dict!=null && dict.getSequence(tokens[0])==null)
						{
						warning("undefined chromosome  in "+bedFile+" "+line);
						continue;
						}
					intervals.put(tokens[0],Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]), true);
					}
				CloserUtil.close(r);
				}
			
			
			

			int prev_tid=-1;
			List<SuspectRgn> buffer=new ArrayList<SuspectRgn>();
			for(;;)
				{
				
				SAMRecord rec=null;
				
				if(merginIter.hasNext())
					{
					rec=merginIter.next();
					progress.watch(rec);
					if(!accept(rec)) continue;
					if(!isClipped(rec)) continue;
					if(intervals!=null && !intervals.containsOverlapping(rec.getReferenceName(),rec.getAlignmentStart(),rec.getAlignmentEnd())) continue;
					}
				
				if(rec==null || prev_tid!=rec.getReferenceIndex())
					{
					if(!buffer.isEmpty())
						{
						info("Analysing buffer buffer.size: "+buffer.size());
						for(SuspectRgn rgn:buffer)
							{
							/* shall we process it ? check at least one depth is OK 
							boolean process_it=false;
							
							for(Input input:inputs)
								{
								int count[]=rgn.sample2count.get(input.sampleName);
								if(count==null) continue; 
								if((count[0]+count[1])*Math.max(1f-min_fraction_of_clipped_records,min_fraction_of_clipped_records) < min_depth) continue;
								process_it=true;
								}
							
							if(!process_it) continue;*/
							Set<Allele> alleles=new HashSet<Allele>();
							int whole_count[]=new int[]{0,0};
							int sum_depth=0;
							List<Genotype> genotypes=new ArrayList<Genotype>();
							int max_depth_here=0;
							
							for(Input input:inputs)
								{
								int count[]=rgn.sample2count.get(input.sampleName);
								if(count==null) count=new int[]{0,0};
								whole_count[0]+=count[0];
								whole_count[1]+=count[1];
								
								
								int depth=0;
								//get depth at position
								SAMRecordIterator it2= input.samFileReaderDepth.queryOverlapping(
										dict.getSequence(rgn.tid).getSequenceName(),
										rgn.pos,
										rgn.pos+(rounding_pos)
										);
								while(it2.hasNext())
									{
									SAMRecord rec2=it2.next();
									if(!accept(rec2)) continue;
									++depth;
									}
								it2.close();
								
								if(depth==0) continue;
								sum_depth+=depth;
								max_depth_here=Math.max(max_depth_here, depth);
								
								float fraction_of_clipped=(float)(count[0]+count[1])/(float)depth;
								info(input.sampleName+" "+rgn+" depth="+(count[0]+count[1])+"/"+depth);
								GenotypeBuilder gb=new GenotypeBuilder(input.sampleName);
								gb.DP(count[0]+count[1]);
								if(depth>=min_depth)
									{
									if(fraction_of_clipped<=0.25)
										{
										gb.alleles(Arrays.asList(reference_allele,reference_allele));
										alleles.add(reference_allele);
										}
									else if(fraction_of_clipped>0.75)
										{
										gb.alleles(Arrays.asList(alternate_allele,alternate_allele));
										alleles.add(alternate_allele);
										}
									else
										{
										gb.alleles(Arrays.asList(reference_allele,alternate_allele));
										alleles.add(reference_allele);
										alleles.add(alternate_allele);
										}
									
									gb.attribute("CN", count[0]+"/"+count[1]);
									}
								else
									{
									//unavailable genotype
									}
								genotypes.add(gb.make());
								}
							if(genotypes.isEmpty())
								{
								info("No genotype "+rgn);
								continue;
								}
							if(max_depth_here< this.min_depth) 
								{
								info("Low max-depth "+max_depth_here+" "+rgn);
								continue;
								}
							if(alleles.size()<=1)
								{
								info("All homozygotes "+rgn);
								continue;//all homozygotes
								}
							Map<String,Object> atts=new HashMap<String,Object>();
							VariantContextBuilder vcb=new VariantContextBuilder();
							vcb.chr(dict.getSequence(rgn.tid).getSequenceName());
							vcb.start(rgn.pos);
							vcb.stop(rgn.pos+this.rounding_pos);
							vcb.alleles(alleles);
							atts.put("DP", sum_depth);
							
							vcb.attributes(atts);
							vcb.genotypes(genotypes);
							
							
							}
						
						buffer.clear();
						}
					if(rec==null)
						{
						break;
						}
					prev_tid=rec.getReferenceIndex();
					}
				
				//this read must be soft clipped in 5' or 3'
				
				
				String sampleName=rec.getReadGroup().getSample();
				
			
				for(int side=0;side<2;++side)
					{
					if(!isClipped(rec,side))continue;
					String clippedSeq=getClippedSequence(rec,side);
					
					
					if(ignore_poly_x)
						{
						int y=1;
						while(y< clippedSeq.length() &&
								same(clippedSeq.charAt(0),clippedSeq.charAt(y)))
							{
							++y;
							}
						if(y==clippedSeq.length() ) continue;
						
						Counter<String> dinucl=new Counter<String>();
						for(y=0;y+1< clippedSeq.length();++y)
							{
							dinucl.incr(clippedSeq.substring(y, y+2));
							}
						if(dinucl.getCountCategories()<=2) continue;
						}					
					int pos=(side==0?rec.getAlignmentStart():rec.getAlignmentEnd());
					
					pos= pos - pos%rounding_pos;
					boolean found_in_buffer=false;
					for(int x=buffer.size()-1;x>=0;--x )
						{
						SuspectRgn rgn=buffer.get(x);
						if(rgn.pos < rec.getUnclippedStart()) break;
						if(rgn.tid!=rec.getReferenceIndex()) continue;
						if(rgn.pos!=pos) continue;
						if(!rgn.sample2count.containsKey(sampleName))
							{
							rgn.sample2count.put(sampleName,new int[]{0,0});
							}
						rgn.sample2count.get(sampleName)[side]++;
						found_in_buffer=true;
						break;
						}
					
					if(!found_in_buffer)
						{
						SuspectRgn rgn=new SuspectRgn();
						rgn.tid=rec.getReferenceIndex();
						rgn.pos=pos;
						rgn.sample2count.put(sampleName,new int[]{(side==0?1:0),(side==0?0:1)});
						buffer.add(rgn);
						if(buffer.size()%10000==0)
							{
							info("buffer size "+buffer.size());
							}
						}
					}

				}
			
			merginIter.close();
			progress.finish();
			return 0;
			}

		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			for(Input input:inputs)
				{
				CloserUtil.close(input);
				}
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new SamFindClippedRegions().instanceMainWithExit(args);
		}
	}
