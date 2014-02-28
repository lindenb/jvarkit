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
	private int min_clip_length=20;
	private boolean ignore_poly_x=false;
	//private int rounding_pos=5;
	//private float min_fraction_of_clipped_records=0.3f;
	private int min_depth=10;
	
	private static class Count
		{
		int count_clipped[]=new int[]{0,0};
		int count_no_clip=0;
		}
	
	private static class SuspectRgn
		{
		int tid=-1;
		int pos=0;
		//boolean found[]=null;
		
		SuspectRgn()
			{
			}
		@Override
		public String toString() {
			return "tid="+tid+":"+pos;
			}
		}
	
	private class Input
		implements Closeable
		{
		int index;
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
	

	private static boolean closeTo(int pos1,int pos2, int max)
		{
		return Math.abs(pos2-pos1)<=max;
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
		if(rec.getReadFailsVendorQualityCheckFlag()) return false;
		Cigar cigar=rec.getCigar();
		if(cigar==null || cigar.isEmpty()) return false;
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
		final int CLOSE_DISTANCE=5;
		int readLength=150;
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
				input.index=inputs.size();
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
			Allele alternate_alleles[]=new Allele[]{ 
				Allele.create("<CLIP5>",false),
				Allele.create("<CLIP3>",false)
				};

			Set<VCFHeaderLine> vcfHeaderLines=new HashSet<VCFHeaderLine>();
			for(Allele alt:alternate_alleles)
				{
				vcfHeaderLines.add(new VCFSimpleHeaderLine(
						"<ID="+alt.getDisplayString()+",Description=\"StructVar\">",
						VCFHeaderVersion.VCF4_1,
						VCFConstants.ALT_HEADER_START.substring(2),
						Arrays.asList("ID", "Description")
						));
				}
			vcfHeaderLines.add(new VCFInfoHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Approximate read depth."));
			vcfHeaderLines.add(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			vcfHeaderLines.add(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			for(int side=0;side<2;++side)
				{
				vcfHeaderLines.add(new VCFFormatHeaderLine("CN"+(side==0?5:3),1,VCFHeaderLineType.Integer,"count clipped in "+(side==0?5:3)+"'"));
				}
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
			
			
			
			CigarElement cigarElements[]=new CigarElement[]{null,null};
			int prev_tid=-1;
			List<SuspectRgn> buffer=new ArrayList<SuspectRgn>();
			for(;;)
				{
				
				SAMRecord rec=null;
				Cigar cigar=null;
				if(merginIter.hasNext())
					{
					rec=merginIter.next();
					progress.watch(rec);
					if(!accept(rec)) continue;
					cigar=rec.getCigar();
					if(cigar==null || cigar.numCigarElements()<2) continue;
					cigarElements[0]=cigar.getCigarElement(0);
					cigarElements[1]=cigar.getCigarElement(cigar.numCigarElements()-1);
					//read must be clipped on 5' or 3' with a good length
					if(!(
						(cigarElements[0].getOperator().equals(CigarOperator.S) && cigarElements[0].getLength()>= this.min_clip_length ) ||
						(cigarElements[1].getOperator().equals(CigarOperator.S) && cigarElements[1].getLength()>= this.min_clip_length )
						))
						{
						continue;
						}
					//not included in bed file ?
					if(intervals!=null && !intervals.containsOverlapping(rec.getReferenceName(),rec.getAlignmentStart(),rec.getAlignmentEnd())) continue;
					}
				
				if(rec==null || prev_tid!=rec.getReferenceIndex())
					{
					if(!buffer.isEmpty())
						{
						info("Analysing buffer buffer.size: "+buffer.size()+". tid:"+prev_tid);
						for(int rgn_idx=0; rgn_idx< buffer.size(); ++rgn_idx )
							{
							if(rgn_idx%100==0) info(""+rgn_idx+"/"+buffer.size());
							SuspectRgn rgn=buffer.get(rgn_idx);
							/* collect the depth of this SuspectRgn for all samples */
							List<Count> counts=new ArrayList<Count>(inputs.size());
							for(Input input:inputs)
								{
								Count count=new Count();
								counts.add(count);
								//if(!rgn.found[input.index]) continue;
								
								
								//get depth at position
								SAMRecordIterator it2= input.samFileReaderDepth.queryOverlapping(
										dict.getSequence(rgn.tid).getSequenceName(),
										Math.max(1, rgn.pos-readLength),
										rgn.pos+readLength
										);
								while(it2.hasNext())
									{
									SAMRecord rec2=it2.next();
									if(!accept(rec2)) continue;
									/* 5' */
									if( rec2.getUnclippedStart()<=rgn.pos &&
										rec2.getUnclippedStart()< rec2.getAlignmentStart() && //read is clipped in 5' 
										closeTo(rec2.getAlignmentStart(), rgn.pos,CLOSE_DISTANCE) &&
										rec2.getAlignmentEnd()>rgn.pos)
										{
										/* this read is clipped in 5' */
										count.count_clipped[0]++;
										}
									/* 3' */
									else if( rec2.getUnclippedEnd()>=rgn.pos &&
											 rec2.getUnclippedEnd()> rec2.getAlignmentEnd() && //read is clipped in 3' 
											closeTo(rec2.getAlignmentEnd(), rgn.pos,CLOSE_DISTANCE) &&
											rec2.getAlignmentStart()<rgn.pos)
										{
										/* this read is clipped in 5' */
										count.count_clipped[1]++;
										}
									else if(rec2.getAlignmentStart()<=rgn.pos &&
											rec2.getAlignmentEnd()>rgn.pos)
										{
										count.count_no_clip++;
										}
									}
								it2.close();
								}/* end foreach input */
							
							int sum_depth=0;
							int max_depth_here=0;
							List<Genotype> genotypes=new ArrayList<Genotype>();	
							Set<Allele> all_alleles=new HashSet<Allele>();
							all_alleles.add(reference_allele);
							boolean found_one_alt=false;
							for(Input input:inputs)
								{
								Count count=counts.get(input.index);
								Set<Allele> sample_alleles=new HashSet<Allele>(3);
								GenotypeBuilder gb=new GenotypeBuilder(input.sampleName);
								gb.DP(count.count_no_clip+count.count_clipped[0]+count.count_clipped[1]);

								for(int side=0;side<2;++side)
									{
									int depth= count.count_no_clip+count.count_clipped[side];
									if(depth==0) continue;//prevent div by 0
									max_depth_here=Math.max(max_depth_here, depth);
									sum_depth+=depth;
									float fraction_of_clipped=(float)(count.count_clipped[side])/(float)depth;
									//info(input.sampleName+" "+rgn+" depth="+(count.count_clipped[side])+"/"+depth);
								
									if(depth>=min_depth)
										{
										if(fraction_of_clipped<=0.25)
											{
											sample_alleles.add(reference_allele);
											}
										else if(fraction_of_clipped>0.75)
											{
											sample_alleles.add(alternate_alleles[side]);
											found_one_alt=true;
											}
										else
											{
											sample_alleles.add(reference_allele);
											sample_alleles.add(alternate_alleles[side]);
											found_one_alt=true;
											}
										
										gb.attribute("CN"+(side==0?5:3),count.count_clipped[side]);
										
										}
									}//end loop side
								if(sample_alleles.isEmpty()) continue;
								gb.alleles(new ArrayList<Allele>(sample_alleles));
								all_alleles.addAll(sample_alleles);
								genotypes.add(gb.make());
								}
							
							if(genotypes.isEmpty())
								{
								//info("No genotype "+rgn);
								continue;
								}
							if(max_depth_here< this.min_depth) 
								{
								//info("Low max-depth "+max_depth_here+" "+rgn);
								continue;
								}
							if(!found_one_alt)
								{
								//info("All homozygotes "+rgn);
								continue;//all homozygotes
								}
							Map<String,Object> atts=new HashMap<String,Object>();
							VariantContextBuilder vcb=new VariantContextBuilder();
							vcb.chr(dict.getSequence(rgn.tid).getSequenceName());
							vcb.start(rgn.pos);
							vcb.stop(rgn.pos);
							vcb.alleles(all_alleles);
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
				
				
				
			
				for(int side=0;side<2;++side)
					{
					CigarElement ce=cigarElements[side];
					if(!ce.getOperator().equals(CigarOperator.S)) continue;
					if(ce.getLength() < this.min_clip_length) continue; 
					String clippedSeq;
					if(side==0)
						{
						clippedSeq=rec.getReadString().substring(0,ce.getLength());
						}
					else
						{
						clippedSeq=rec.getReadString().substring(rec.getReadLength()-ce.getLength());
						}
					
					
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
					
					boolean found_in_buffer=false;
					for(int x=buffer.size()-1;x>=0;--x )
						{
						SuspectRgn rgn=buffer.get(x);
						if(rgn.pos < rec.getUnclippedStart() - readLength ) break;
						if(rgn.tid!=rec.getReferenceIndex()) continue;
						if(rgn.pos!=pos) continue;
						//rgn.found[sample2input.get(sampleName).index]=true;
						found_in_buffer=true;
						break;
						}
					
					if(!found_in_buffer)
						{
						SuspectRgn rgn=new SuspectRgn();
						rgn.tid=rec.getReferenceIndex();
						rgn.pos=pos;
						//rgn.found=new boolean[inputs.size()];
						//Arrays.fill(rgn.found, false);
						//rgn.found[sample2input.get(sampleName).index]=true;
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
