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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimplePosition;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;

/**
 BEGIN_DOC
 
## Input

The input is a file containing a list of path to the bam.
 
## Example

```
$ find ./testdata/ -type f -name "*.bam" | \
 java -jar dist/findallcoverageatposition.jar -p rotavirus:100


#File              CHROM      POS  SAMPLE  DEPTH  M    I  D  N  S   H  P  EQ  X  Base(A)  Base(C)  Base(G)  Base(T)  Base(N)  Base(^)  Base(-)
./testdata/S4.bam  rotavirus  100  S4      126    126  0  0  0  29  0  0  0   0  5        0        0        121      0        0        0
./testdata/S1.bam  rotavirus  100  S1      317    317  1  0  0  50  0  0  0   0  27       0        1        289      0        1        0
./testdata/S2.bam  rotavirus  100  S2      311    311  0  1  0  60  0  0  0   0  29       1        0        281      0        0        1
./testdata/S3.bam  rotavirus  100  S3      446    446  1  0  0  86  0  0  0   0  39       0        1        406      0        1        0

```

## See also

 * [https://twitter.com/pjacock/status/538300664334798848](https://twitter.com/pjacock/status/538300664334798848)
 * [https://twitter.com/yokofakun/status/538300434109456385](https://twitter.com/yokofakun/status/538300434109456385)
 * [https://twitter.com/pjacock/status/538299549455233024](https://twitter.com/pjacock/status/538299549455233024)
 * FindAVariation

## History

 + 20210303 : bug fixed in https://github.com/lindenb/jvarkit/issues/180


END_DOC
 */
@Program(name="findallcoverageatposition",
	keywords={"bam","coverage","search","depth"},
	description="Find depth at specific position in a list of BAM files. My colleague Estelle asked: in all the BAM we sequenced, can you give me the depth at a given position ?",
	biostars= {259223,250099,409942},
	modificationDate="20210818",
	creationDate="20141128"
	)
public class FindAllCoverageAtPosition extends Launcher
	{
	private static final Logger LOG = Logger.build(FindAllCoverageAtPosition.class).make();
	private static final char INSERTION_CHAR='^';
	private static final char DELETION_CHAR='-';
	private static final char BASES_To_PRINT[]=new char[]{'A','C','G','T','N',INSERTION_CHAR,DELETION_CHAR};

	@Parameter(names={"-p","--position"},description="-p chrom:pos . Multiple separated by space. Add this chrom/position. Required")
	private List<String> positionStrs = new ArrayList<String>();
	@Parameter(names={"-f","--posfile"},description="File containing positions. if file suffix is '.bed': all positions in the range will be scanned.")
	private Path positionFile = null;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--groupby","--partition"},description="Group Reads by. "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition groupBy=SAMRecordPartition.sample;
	@Parameter(names={"-filter","--filter"},description="[20171201](moved to jexl). "+SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamRecordJEXLFilter.StringConverter.class)
	private SamRecordFilter filter = SamRecordJEXLFilter.buildDefault();
	@Parameter(names={"-r","-R","--reference"},description="[20171201]"+Launcher.INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path referenceFileFile=null;
	@Parameter(names={"-x","--extend"},description="[20190218]extend by 'x' base to try to catch close with clipped reads. " + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int extend=500;
	@Parameter(names={"-Q","--mapq"},description="Min mapping quality. Dicard reads having MAPQ < 'x'")
	private int mapq=1;
	@Parameter(names={"-clip","--clip"},description="use clipped bases.")
	private boolean use_clipped_bases = false;

	
	
	private ReferenceSequenceFile indexedFastaSequenceFile=null;
	private GenomicSequence genomicSequence=null;
	private SAMSequenceDictionary fastaDict = null;
	
	private static class CigarAndBases
		{
		Counter<CigarOperator> operators = new Counter<>();
		Counter<Character> bases = new Counter<>();
		}
	private PrintWriter out=null;
	private SamReaderFactory samReaderFactory;
    public FindAllCoverageAtPosition()
    	{
    	}		
    
    private SimplePosition convertFromSamHeader(final Path f,SAMFileHeader h,final SimplePosition src)
    	{
    	final SAMSequenceDictionary dict=h.getSequenceDictionary();
    	if(dict==null)
    		{
    		LOG.warn("No dictionary in "+f);
    		return null;
    		}
    	final ContigNameConverter converter = ContigNameConverter.fromOneDictionary(dict); 
    	final String ctg =  converter.apply(src.getContig());
    	if(ctg==null) return null;
    	return src.renameContig(ctg);
    	}
    
    private char getReferenceAt(final String contig,int pos1) {
    	if(this.indexedFastaSequenceFile==null) return '.';
    	if(this.genomicSequence==null || !this.genomicSequence.getChrom().equals(contig)) {
    		final SAMSequenceDictionary dict= this.indexedFastaSequenceFile.getSequenceDictionary();
    		if(dict==null) return '.';
    		
    		final ContigNameConverter converter = ContigNameConverter.fromOneDictionary(dict); 
        	final String newctg =  converter.apply(contig);
        	final SAMSequenceRecord rec= (newctg==null?null:dict.getSequence(newctg));
        	if(rec!=null) genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, newctg);
    		}
    	if(genomicSequence==null ||pos1<0 || pos1>genomicSequence.length()) return '.';
    	return genomicSequence.charAt(pos1-1);
    	}

    private void scan(final BufferedReader in,final Set<SimplePosition> mutations) throws Exception
    	{

    	final String DEFAULT_SAMPLE_NAME="(undefined)";
    	String line;
    	while((line=in.readLine())!=null)
			{
    		if(out.checkError()) break;
			if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
			final Path f= Paths.get(line);
			if(!Files.exists(f)) continue;
			if(!Files.isRegularFile(f)) continue;
			if(!Files.isReadable(f)) continue;
			final String filename=f.getFileName().toString();
			
			if(filename.endsWith(FileExtensions.CRAM)) {
				if(referenceFileFile==null) {
					LOG.warn("skipping "+f+" not reference was provided");
					continue;
				}
				final SAMSequenceDictionary dict  = SAMSequenceDictionaryExtractor.extractDictionary(f);
				if(dict==null) continue;
				if(!SequenceUtil.areSequenceDictionariesEqual(dict,this.fastaDict)) {
					LOG.warn("skipping "+f+" because this is not the same sequence dictionary");
					continue;
				}
			} 
			
			if(!(filename.endsWith(FileExtensions.BAM) || 
				filename.endsWith(FileExtensions.CRAM))) continue;
						
    			
		
			try
				{
				try(SamReader samReader = this.samReaderFactory.open(f)) {
					if(!samReader.hasIndex())
						{
						LOG.warn("no index for "+f);
						continue;
						}
					final SAMFileHeader header=samReader.getFileHeader();
					for(final SimplePosition src:mutations)
						{
						final Map<String, CigarAndBases> sample2count=new TreeMap<String,CigarAndBases>();
						for(final SAMReadGroupRecord rg:header.getReadGroups())
							{
							if(rg!=null)
								{
								String sn=this.groupBy.apply(rg);
								if(sn!=null && !sn.trim().isEmpty())
									{
									sample2count.put(sn, new CigarAndBases());
									}
								}
							}
						
						if(sample2count.isEmpty())
							{
							sample2count.put(DEFAULT_SAMPLE_NAME, new CigarAndBases());
							}
	
						
						
						final SimplePosition m = convertFromSamHeader(f,header,src);
						if(m==null) continue;
						try(CloseableIterator<SAMRecord> iter = samReader.query(m.getContig(),
								Math.max(1,m.getPosition()-this.extend),
								m.getPosition()+this.extend,
								false
								)) {
							while(iter.hasNext())
								{
								final SAMRecord rec=iter.next();
								if(rec.getReadUnmappedFlag()) continue;
								if(rec.getMappingQuality() < this.mapq) continue;
								if(this.filter.filterOut(rec)) continue;
								final Cigar cigar=rec.getCigar();
								if(cigar==null || cigar.isEmpty()) continue;
								if(rec.getUnclippedEnd() < m.getPosition()) continue;
								if(rec.getUnclippedStart() > m.getPosition()) continue;
								if(rec.getReadBases()==SAMRecord.NULL_SEQUENCE) continue;
								
								final String readString = rec.getReadString().toUpperCase();
								final String sampleName= this.groupBy.getPartion(rec,DEFAULT_SAMPLE_NAME);
								
								CigarAndBases counter= sample2count.get(sampleName);
								if(counter==null)
									{
									counter=new CigarAndBases();
									sample2count.put(sampleName, counter);
									}	
								
								
								int ref= rec.getUnclippedStart();
								int readPos = 0;
								for(int k=0;k<cigar.numCigarElements() && ref< m.getPosition()+1;++k)
									{
									final CigarElement ce=cigar.getCigarElement(k);
									final CigarOperator op=ce.getOperator();
									switch(op)
										{
										case P: break;
										case I: 
											{
											if(ref==m.getPosition())
												{
												counter.operators.incr(op);
												counter.bases.incr(INSERTION_CHAR);
												}
											readPos += ce.getLength();
											break;
											}
										case D:case N:
										case M: case X: case EQ: 
										case H: case S:
											{
											for(int i=0;i< ce.getLength();++i )
												{
												if(ref==m.getPosition())
													{
													counter.operators.incr(op);
													switch(op)
														{
														case H:
															{
															if(use_clipped_bases) counter.bases.incr('N');
															break;
															}
														case S:
															{
															if(use_clipped_bases) counter.bases.incr(readString.charAt(readPos));
															break;
															}
														case M:case X:case EQ:
															counter.bases.incr(readString.charAt(readPos));
															break;
														case D:case N:
															counter.bases.incr(DELETION_CHAR);
															break;
														default:break;
														}
													}	
												if(op.consumesReadBases()) ++readPos;
												ref++;
												}
											break;
											}
										default: throw new RuntimeException("unknown operator:"+op);
										}
									}
								}
						}
							
						
						for(final String sample:sample2count.keySet())
							{
							final CigarAndBases counter= sample2count.get(sample);
							
							out.print(f);
							out.print('\t');
							out.print(m.getContig());
							out.print('\t');
							out.print(m.getPosition());
							
							if(this.indexedFastaSequenceFile!=null) {
								out.print('\t');
								out.print(getReferenceAt(m.getContig(),m.getPosition()));
								}
							
							out.print('\t');
							out.print(sample);
							out.print('\t');
							out.print(
									counter.operators.count(CigarOperator.M)+
									counter.operators.count(CigarOperator.EQ)+
									counter.operators.count(CigarOperator.X)+
									(use_clipped_bases?counter.operators.count(CigarOperator.H)+counter.operators.count(CigarOperator.S):0)
									);
							for(final CigarOperator op:CigarOperator.values())
								{
								out.print('\t');
								out.print(counter.operators.count(op));
								}
							for(char c:BASES_To_PRINT)
								{
								out.print('\t');
								out.print(counter.bases.count(c));
								}
							
							out.println();
							}
							
						}//end of loop mutations
				}//end of samReader
				}
			catch(final Throwable err)
				{
				LOG.error(err);
				throw err;
				}
			finally
				{
				}    				
			}		
	    
    	}
    
    @Override
    public int doWork(final List<String> args) {
    	final Set<SimplePosition> mutations=new TreeSet<>();
    	if(this.extend<1) {
    		LOG.error("-x should be >=1");
    		return -1;
    	}
		
		try
			{
			this.samReaderFactory=  SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);

			
			if(this.referenceFileFile!=null) {
				this.indexedFastaSequenceFile = htsjdk.samtools.reference.ReferenceSequenceFileFactory.getReferenceSequenceFile(this.referenceFileFile);
				this.samReaderFactory.referenceSequence(this.referenceFileFile);
				this.fastaDict = SequenceDictionaryUtils.extractRequired(this.indexedFastaSequenceFile);
				}
			
			this.positionStrs.
				stream().
				flatMap(S->Arrays.stream(S.split("[ \t;]+"))).
				filter(S->!StringUtils.isBlank(S)).
				map(S->new SimplePosition(S)).
				forEach(X->mutations.add(X));
			
			
			if(this.positionFile!=null) {
				if(	this.positionFile.toString().endsWith(".bed") || 
					this.positionFile.toString().endsWith(".bed.gz")) {
					try(BedLineReader br=new BedLineReader(this.positionFile)) {
						br.stream().
						filter(B->B!=null).
						forEach(bedLine->{
							for(int x=bedLine.getStart();x<=bedLine.getEnd();++x)
								{
								mutations.add(new SimplePosition(bedLine.getContig(),x));
								}
							});
						}
					}
				else
					{
					try(BufferedReader r2= IOUtils.openPathForBufferedReading(this.positionFile)) {
						r2.lines().
							filter(L->!StringUtils.isBlank(L)).
							filter(L->!L.startsWith("#")).
							forEach(L->mutations.add(new SimplePosition(L)));
						}
					}
				
			}
		
		
			if(mutations.isEmpty())
				{
				LOG.fatal( "undefined position \'str\'");
				return -1;
				}
		
			LOG.info("number of mutations "+mutations.size());
			
			
			this.out=this.openPathOrStdoutAsPrintWriter(this.outputFile);
			
			out.print("#File");
			out.print('\t');
			out.print("CHROM");
			out.print('\t');
			out.print("POS");

			if(this.indexedFastaSequenceFile!=null)
				{
				out.print('\t');
				out.print("REF");
				}
			
			out.print('\t');
			out.print(this.groupBy.name().toUpperCase());
			out.print('\t');
			out.print("DEPTH");
			for(final CigarOperator op:CigarOperator.values())
				{
				out.print('\t');
				out.print(op.name());
				}
			for(char c:BASES_To_PRINT)
				{
				out.print('\t');
				out.print("Base("+c+")");
				}
			
			out.println();

			
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				try(BufferedReader r = new BufferedReader(new InputStreamReader(stdin()))) {
					scan(r,mutations);
					}
				}
			else
				{				
				for(final String filename: args)
					{
					LOG.info("Reading from "+filename);
					try(BufferedReader r=IOUtils.openURIForBufferedReading(filename)) {
						scan(r,mutations);
						}
					}
				}
			this.out.flush();
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.severe(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			CloserUtil.close(this.out);
			}
		}
	
	/**
	 * main
	 */
	public static void main(String[] args) {
		new FindAllCoverageAtPosition().instanceMainWithExit(args);

	}

}
