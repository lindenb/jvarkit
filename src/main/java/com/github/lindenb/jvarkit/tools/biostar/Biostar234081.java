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

package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.CloserUtil;
/**
BEGIN_DOC

## Example

```
 $ cat toy.sam 
@SQ SN:ref  LN:45
@SQ SN:ref2 LN:40
r001    163 ref 7   30  1M2X5=4I4M1D3M  =   37  39  TTAGATAAAGAGGATACTG*XX:B:S,12561,2,20,112

 $ java -jar dist/biostar234081.jar toy.sam 
@HD VN:1.5  SO:unsorted
@SQ SN:ref  LN:45
@SQ SN:ref2 LN:40
r001    163 ref 7   30  8M4I4M1D3M  =   37  39  TTAGATAAAGAGGATACTG*XX:B:S,12561,2,20,112
```

END_DOC
 */
@Program(name="biostar234081",
	description="convert extended CIGAR to regular CIGAR ('X','=' -> 'M')",
	keywords={"sam","bam","cigar"},
	biostars=234081
	)
public class Biostar234081 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar234081.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
	@Override
	public int doWork(final List<String> args) { 
	SamReader in =null;
	SAMFileWriter w=null;  
	SAMRecordIterator iter=null;  
	try {
		in = super.openSamReader(oneFileOrNull(args));
		final SAMFileHeader header= in.getFileHeader();
		JVarkitVersion.getInstance().addMetaData(this, header);
		final SAMProgramRecord prg = header.createProgramRecord();
		prg.setProgramName(this.getProgramName());
		prg.setProgramVersion(this.getGitHash());
		prg.setCommandLine(this.getProgramCommandLine());
		w = this.writingBamArgs.openSAMFileWriter(this.outputFile,header,true);
		final ProgressFactory.Watcher<SAMRecord> progress = ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
		iter=in.iterator();
		while(iter.hasNext())
			{
			final SAMRecord rec = progress.apply(iter.next());
			if(!rec.getReadUnmappedFlag() &&
					rec.getCigar()!=null &&
					rec.getCigar().getCigarElements().
						stream().
						map(C->C.getOperator()).
						anyMatch(OP->(OP.equals(CigarOperator.EQ) || OP.equals(CigarOperator.X)))
						)
				{
				final Cigar cigar = rec.getCigar();
				final List<CigarElement> elements = new ArrayList<>(cigar.numCigarElements());

				for(int i=0;i< cigar.numCigarElements();++i)
					{
					CigarElement ce =cigar.getCigarElement(i);
					switch(ce.getOperator())
						{
						case X:
						case EQ: ce =new CigarElement(ce.getLength(), CigarOperator.M);break;
						default:break;
						}
					if(!elements.isEmpty() && elements.get(elements.size()-1).getOperator()==ce.getOperator())
						{
						elements.set(elements.size()-1,
								new CigarElement(
										ce.getLength()+elements.get(elements.size()-1).getLength(),
										ce.getOperator())
								);
						}
					else
						{
						elements.add(ce);
						}
					}
				rec.setCigar(new Cigar(elements));
				rec.setAttribute(SAMTag.PG.name(), prg.getId());
				}
			w.addAlignment(rec);
			}
		progress.close();
		iter.close();iter=null;
		w.close();w=null;
		return RETURN_OK;
	} catch (final Exception err) {
		LOG.error(err);
		return -1;
	}	
	finally {
		CloserUtil.close(iter);
		CloserUtil.close(in);
		CloserUtil.close(w);
		}
	}

	public static void main(final String[] args)
		{
		new Biostar234081().instanceMainWithExit(args);
		}

	}
