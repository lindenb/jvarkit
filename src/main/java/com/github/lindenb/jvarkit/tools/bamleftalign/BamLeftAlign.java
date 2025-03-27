package com.github.lindenb.jvarkit.tools.bamleftalign;

import java.util.Collections;
import java.util.List;
import java.util.function.Function;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassBamLauncher;
import com.github.lindenb.jvarkit.samtools.SAMRecordLeftAligner;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;


/**
BEGIN_DOC

## Example

```

```


END_DOC
 */
@Program(name="bamleftalign",
	keywords={"bam","deletion","cram"},
	description="Left Align Reads around deletions",
	modificationDate="20250327",
	creationDate="20250327",
	jvarkit_amalgamion = true,
	menu="BAM Manipulation"
	)
public class BamLeftAlign extends OnePassBamLauncher {
	private static final Logger LOG = Logger.build(BamLeftAlign.class).make();
	private enum What {none,only,discard};
	@Parameter(names={"--filter"},description="none: keep any read (realigned or not); only: only keep realigned reads; discard: discard realigned reads.")
	protected  What what_to_do = What.none;

	
	private ReferenceSequenceFile reference = null;
	private SAMRecordLeftAligner leftAligner=null;
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	@Override
	protected int beforeSam() {
		try {
			reference=ReferenceSequenceFileFactory.getReferenceSequenceFile(super.getRequiredReferencePath());
			}
		catch(Throwable err) {
			getLogger().error(err);
			return -1;
			}
		this.leftAligner = new SAMRecordLeftAligner(reference);
		return super.beforeSam();
		}
	
	private List<SAMRecord> leftAlign(SAMRecord rec) {
		rec = this.leftAligner.apply(rec);
		switch(what_to_do) {
			case none: break;
			case discard: if(this.leftAligner.previousRecordWasRealigned()) return Collections.emptyList();break;
			case only: if(!this.leftAligner.previousRecordWasRealigned()) return Collections.emptyList();break;
			default: throw new IllegalStateException();
			}
		return Collections.singletonList(rec);
		}
	
	
	@Override
	protected Function<SAMRecord, List<SAMRecord>> createSAMRecordFunction() {
			return REC->leftAlign(REC);
			}
	
	@Override
	protected void afterSam() {
		if(leftAligner!=null) {
			leftAligner.dispose();
			leftAligner=null;
			}
		if(reference!=null) {
			try { reference.close();}
			catch(Throwable err) {}
			reference=null;
			}
		super.afterSam();
		}
	public static void main(String[] args) {
		new BamLeftAlign().instanceMainWithExit(args);
	}

}
