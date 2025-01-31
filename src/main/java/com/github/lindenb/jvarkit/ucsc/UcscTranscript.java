/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.ucsc;

import java.util.AbstractList;
import java.util.Collections;
import java.util.List;
import java.util.OptionalInt;
import java.util.stream.Collectors;

import org.apache.commons.jexl2.JexlContext;

import com.github.lindenb.jvarkit.bed.BedInterval;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.util.bio.KozakSequence;

import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;

/**
 * Describe a UCSC transcript in the genepred format
 * @author lindenb
 *
 */
public interface UcscTranscript extends Feature , BedInterval{
/** get Transcript Id */
public String getTranscriptId();
/** get name2 or 'gene' if available */
public String getName2();
public Strand getStrand();
public boolean isPositiveStrand();
public boolean isNegativeStrand();
public int getTxStart();
public int getTxEnd();

/** return CDS start or throw an exception if it is not a protein coding transcript */
public int getCdsStart();
/** return CDS end or throw an exception if it is not a protein coding transcript */
public int getCdsEnd();


/** return sum of length on reference of exons */
public default int getTranscriptLength() {
	return getExons().stream().mapToInt(E->E.getLengthOnReference()).sum();
	}

/** get TxStart */
@Override
public default int getStart() {
	return getBedStart() +1 ;
	}

/** get TxEnd */
@Override
public default  int getEnd() {
	return getBedEnd();
	}
/** get count of exons */
public int getExonCount();

public default boolean hasIntrons() {
	return this.getExonCount() > 1;
}

public default int getIntronCount() {
	return Math.max(0, getExonCount()-1);
	}

public int getExonStart(int idx);
public int getExonEnd(int idx);


public boolean isProteinCoding();


/** return 0-based genomic position of the **A**TG initiator */
public default int getGenome0ATG() {
	if(!isProteinCoding()) throw new IllegalStateException("non coding protein");
	if(isPositiveStrand()) {
		return getCdsStart()-1;//0 based ATG
		}
	else
		{
		return getCdsEnd()-1;//ATG is *before* the END position
		}
	}


/** negate isCoding */
public default boolean isNonCoding() {
	return !isProteinCoding();
	}
/** alias of isProteinCoding */
public default boolean isCoding() {
	return isProteinCoding();
	}


/** get i-th intron */
public default Intron getIntron(int idx) {
	return new Intron(this,idx);
	}

public default Exon getExon(int idx) {
	return new Exon(this,idx);
	}

public default List<Exon> getExons() {
	return new AbstractList<Exon>() {
		@Override
		public Exon get(int idx) {
			return getExon(idx);
			}
		@Override
		public int size() {
			return getExonCount();
			}
		};
	}



public default List<Intron> getIntrons() {
	if(!hasIntrons()) return Collections.emptyList();
	return new AbstractList<Intron>() {
		@Override
		public Intron get(int idx) {
			return getIntron(idx);
			}
		@Override
		public int size() {
			return getIntronCount();
			}
		};
	}
/** alias for getCDSs() , in the genomic order */
public default List<CDS> getCDS() {
	return getCDSs();
	}

/** return list of CDS, in the genomic order */
public default List<CDS> getCDSs() {
	if(!isProteinCoding()) Collections.emptyList();
	final List<CDS> L = getExons().stream().
		filter(EX->CoordMath.overlaps(
				EX.getStart(), EX.getEnd(),
				getCdsStart(), getCdsEnd()
				)).
		map(EX->new CDS(EX)).
		collect(Collectors.toList());
	return L;
	}




public abstract class Component implements Locatable, BedInterval {
	public abstract UcscTranscript getTranscript();
	
	@Override
	public int hashCode() {
		int i = getContig().hashCode();
		i= i*31 + getTranscript().getTranscriptId().hashCode();
		i= i*31 + Integer.hashCode(getStart());
		i= i*31 + Integer.hashCode(getEnd());
		return i;
		}
	
	
	public boolean containsGenomicLoc1(int pos1) {
		return getStart() <= pos1 && pos1 <= getEnd();
	}
	
	public abstract String getName();
	
	@Override
	public final int getBedStart() {
		return getStart()-1;
		}
	
	@Override
	public final int getBedEnd() {
		return getEnd();
		}
	
	public Strand getStrand() {
		return getTranscript().getStrand();
		}
	
	public boolean isPositiveStrand() {
		return getTranscript().isPositiveStrand();
		}
	
	public boolean isNegativeStrand() {
		return getTranscript().isNegativeStrand();
		}
	
	@Override
	public String getContig() {
		return getTranscript().getContig();
		}
	@Override
	public String toString() {
		return getClass().getSimpleName()+":"+getTranscript().getTranscriptId()+":"+getContig()+":"+getStart()+"-"+getEnd();
		}
	
	public Interval toInterval() {
		return new Interval(getContig(), getStart(), getEnd(), isNegativeStrand(), getName());
		}
	}



/************************************************************************
 * 
 * Exon
 *
 */
public class Exon extends Component {
	private final UcscTranscript owner;
	private final int exon_index;
	Exon(final UcscTranscript owner,int exon_index) {
		this.owner = owner;
		this.exon_index = exon_index;
		}
	
	
	/** return 0-based exon index in the genome order */
	public int getIndex() {
		return this.exon_index;
		}
	
	@Override
	public UcscTranscript getTranscript() {
		return this.owner;
		}
	
	@Override
	public boolean equals(final Object obj) {
		if(this==obj) return true;
		if(obj==null || !(obj instanceof Exon)) return false;
		final Exon o = Exon.class.cast(obj);
		return exon_index == o.exon_index &&
				contigsMatch(o) &&
				getTranscript().getTranscriptId().equals(o.getTranscript().getTranscriptId());
		}

	/** return 1-based exon index , according to strand */
	public int getHumanIndex() {
		return isPositiveStrand()?exon_index+1:getTranscript().getExonCount()-exon_index;
		}
	@Override
	public int getStart() {
		return getTranscript().getExonStart(this.exon_index);
		}
	@Override
	public int getEnd() {
		return getTranscript().getExonEnd(this.exon_index);
		}
	@Override
	public String getName() {
		return getTranscript().getTranscriptId()+".Exon"+ getHumanIndex();
		}
	}



/************************************************************************
 *
 * Intron
 *
 */

public default int getIntronStart(int idx) {
	return getExonEnd(idx) +1;
	}
public default int getIntronEnd(int idx) {
	return getExonStart(idx+1) -1;
	}


public class Intron extends Component {
	private final UcscTranscript owner;
	private final int intron_index;
	private Intron(final UcscTranscript owner,int intron_index) {
		this.owner = owner;
		this.intron_index = intron_index;
		}
	/** return 0-based intron index in the genome order */
	public int getIndex() {
		return this.intron_index;
		}
	@Override
	public UcscTranscript getTranscript() {
		return this.owner;
		}

	/** return 1-based intron index , according to strand */
	public int getHumanIndex() {
		return isPositiveStrand()?intron_index+1:getTranscript().getIntronCount()-intron_index;
		}
	
	@Override
	public int getStart() {
		return getTranscript().getIntronStart(this.intron_index);
		}
	@Override
	public int getEnd() {
		return getTranscript().getIntronEnd(this.intron_index);
		}
	@Override
	public String getName() {
		return getTranscript().getTranscriptId()+".Intron"+getHumanIndex();
		}
	
	@Override
	public boolean equals(final Object obj) {
		if(this==obj) return true;
		if(obj==null || !(obj instanceof Intron)) return false;
		final Intron o = Intron.class.cast(obj);
		return intron_index==o.intron_index && 
				contigsMatch(o) && 
				getTranscript().getTranscriptId().equals(o.getTranscript().getTranscriptId());
		}
	/** get Surrounding exons */
	public Exon[] getEnclosingExons() {
		return new Exon[] {
				getTranscript().getExon(this.intron_index),
				getTranscript().getExon(this.intron_index+1),
			};
		}
	}

/** ===================================================================*/
public abstract class ExonComponent extends Component {
	private final Exon exon;
	private ExonComponent(final Exon exon) {
		this.exon = exon;
		}
	
	
	@Override
	public UcscTranscript getTranscript() {
		return this.exon.getTranscript();
		}
	
	public final Exon getExon() {
		return this.exon;
		}
	@Override
	public String getName() {
		return getTranscript().getTranscriptId()+"."+getClass().getSimpleName() + "."+getStart()+"-"+getEnd();
		}
	}



public class CDS extends ExonComponent {
	private CDS(final Exon exon) {
		super(exon);
		}
	@Override
	public int getStart() {
		return Math.max(
			getExon().getStart(),
			getTranscript().getCdsStart()
			);
		}
	@Override
	public int getEnd() {
		return Math.min(
			getExon().getEnd(),
			getTranscript().getCdsEnd()
			);
		}
	
	
	@Override
	public boolean equals(Object obj) {
		if(this==obj) return true;
		if(obj==null || !(obj instanceof CDS)) return false;
		final CDS o = CDS.class.cast(obj);
		return getStart()==o.getStart() && 
				getEnd()==o.getEnd() && 
				contigsMatch(o) && 
				getTranscript().getTranscriptId().equals(o.getTranscript().getTranscriptId());
		}
	}

/** ===================================================================*/
public abstract class UTR extends ExonComponent {
	protected UTR(final Exon exon) {
		super(exon);
		}
	public abstract boolean isUTR5();
	public abstract boolean isUTR3();
	}

/** return true if mRNA has UTR */
public default boolean hasUTRs() {
	return hasUTR5() || hasUTR3();
}
/** return  list of mixed UTR5/UTR3 */
public List<UTR> getUTRs();

/** return true if there is a 5' UTR */
public default boolean hasUTR5() {
	return isProteinCoding() &&
		(isPositiveStrand()?
			getTxStart() < getCdsStart():
			getCdsEnd() < getTxEnd())
			;
	}


/** return true if there is a 3' UTR */
public default boolean hasUTR3() {
	return isProteinCoding() &&
		(isPositiveStrand()?
			getCdsEnd() < getTxEnd():
			getTxStart() < getCdsStart()
			);
	}


/** ===================================================================*/
public abstract class Codon extends Component {
	private final CDS cds;
	protected Codon(final CDS cds) {
		this.cds=cds;
		}
	public CDS getCDS() { return this.cds;}
	public Exon getExon() { return this.getCDS().getExon();}
	@Override public UcscTranscript getTranscript() { return getExon().getTranscript();}
	public abstract boolean isStartCodon();
	public final boolean isStopCodon() { return !isStartCodon();}
	}


public List<Codon> getCodons();

public interface RNA extends CharSequence , Locatable, BedInterval {
	/** convert 0-based rna to 0-based genomic */
	public int convertToGenomic0Coordinate0(int rnaPos0);
	/** convert the genomic position to the position in the RNA, return empty( if RNA not in genomic pos */
	public OptionalInt convertGenomic0ToRnaCoordinate0(int genomicPos0);
	/** return associated transcript, it it's a micro-ORF, the new transcript associated to the ORF will be returned */
	public UcscTranscript getTranscript();
	public default Strand getStrand() { return getTranscript().getStrand();}
	public default String getContig() { return getTranscript().getContig(); }
	public default boolean isPositiveStrand() { return getTranscript().isPositiveStrand();}
	public default boolean isNegativeStrand() { return getTranscript().isNegativeStrand();}
}


public interface MessengerRNA extends RNA {
	public default boolean hasCodingRNA() {
		return getTranscript().isProteinCoding();
		}
	public default boolean hasUpstreamUntranslatedRNA() {
		return getTranscript().hasUTR5();
		}
	public UntranslatedRNA getUpstreamUntranslatedRNA();
	public default boolean hasDownstreamUntranslatedRNA() {
		return getTranscript().hasUTR3();
		}
	public UntranslatedRNA getDownstreamUntranslatedRNA();
	public CodingRNA getCodingRNA();
	
	}

public interface CodingRNA extends RNA {
	public MessengerRNA getMessengerRNA();
	public KozakSequence getKozakSequence();
	public default KozakSequence.Strength getKozakStrength() {
		return getKozakSequence().getStrength();
		}
	/** convert cDNA position to mRNA 0-based position */
	public int convertCoding0ToMessenger0(int cDNA0);
 	public Peptide getPeptide();
	}


public interface UntranslatedRNA extends RNA {
	public MessengerRNA getMessengerRNA();
	/* find micro ORF starting in this UTR */
	public List<CodingRNA> getORFs();
	}



/** return  mRNA for this chromosome */
public MessengerRNA getMessengerRNA(final CharSequence chromosomeSequence);

/** return  mRNA for this chromosome with a fake genomic sequence. base of the genome should never be asked*/
public default MessengerRNA getMessengerRNA() {
	return getMessengerRNA(new AbstractCharSequence() {	
		@Override
		public int length() {
			return Integer.MAX_VALUE-100;
		}
		@Override
		public char charAt(int arg0) {
			return 'N';
		}
	});
}



public interface Peptide extends CharSequence {
	public CodingRNA getCodingRNA();
	public int[] convertToGenomic0Coordinates(int pepPos0);
	/** convert the genomic position to the position in the peptide, return -1 if peptide not in genomic pos */
	public default int convertGenomic0ToPeptideCoordinate0(int genomicPos0)
		{
		final int rnaIdx=getCodingRNA().convertToGenomic0Coordinate0(genomicPos0);
		if(rnaIdx==-1) return  -1;
		return rnaIdx/3;
		}
	}



/** return this transcript as a JEXL context */
public default JexlContext asJexlContext() {
	return new AsJexlContext(this);
	}

static class AsJexlContext implements JexlContext
	{
	private final UcscTranscript kg;
	AsJexlContext(final UcscTranscript kg) {
		this.kg=kg;
		}
	@Override
	public Object get(final String s) {
		if(s.equals("kg") || s.equals("gene") || s.equals("transcript")) return this.kg;
		return null;
		}
	@Override
	public boolean has(final String s) {
		if(s.equals("kg") || s.equals("gene") || s.equals("transcript")) return true;
		return false;
		}
	@Override
	public void set(final String arg0, final Object arg1) {
		throw new UnsupportedOperationException("cannot set "+arg0+" to "+this.kg);
		}
	@Override
	public String toString() {
		return kg.toString();
		}
	}

}
