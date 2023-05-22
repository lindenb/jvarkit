package com.github.lindenb.jvarkit.ucsc;

import java.util.AbstractList;
import java.util.List;
import java.util.OptionalInt;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.ucsc.UcscTranscript.UTR5;

import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.annotation.Strand;

public interface UcscTranscript extends Locatable {
	public OptionalInt getCdsStart();
	public OptionalInt getCdsEnd();
	public Strand getStrand();
	public int getExonCount();
	public AttributeMap getAttributes();
	public default int getIntronCount() {
		return getExonCount() - 1;
		}
	public UcscTranscript.Exon getExon(int idx);
	public default List<Exon> getExons() {
		return new AbstractList<Exon>()
			{
			@Override
			public int size()
				{
				return getExonCount();
				}
			@Override
			public Exon get(int index)
				{
				return getExon(index);
				}
			};
		}
	
	public UcscTranscript.Intron getIntron(int idx);
	public default List<Intron> getIntrons() {
		return new AbstractList<Intron>()
			{
			@Override
			public int size()
				{
				return getIntronCount();
				}
			@Override
			public Intron get(int index)
				{
				return getIntron(index);
				}
			};
		}
	

	
	public default boolean isPositiveStrand() {
		return getStrand().equals(Strand.FORWARD);
		}
	public default boolean isNegativeStrand() {
		return getStrand().equals(Strand.REVERSE);
		}
	public default boolean hasCDS() {
		return getCdsStart().isPresent() && getCdsEnd().isPresent();
		}

	public List<UTR5> getUTR5AsList();
	public List<UTR3> getUTR3AsList();
	
	public interface TranscriptComponent extends Locatable {
		public UcscTranscript getTranscript();
		@Override
		public default String getContig() { return getTranscript().getContig();}
		public default Strand getStrand()  { return getTranscript().getStrand();}
		public default boolean isPositiveStrand() {
			return getTranscript().isPositiveStrand();
			}
		public default boolean isNegativeStrand() {
			return getTranscript().isNegativeStrand();
			}
		public default IntStream getGenomicPositions() {
			return IntStream.rangeClosed(getStart(), getEnd());
			}
		
		public String getName();
		}
	public interface Exon extends TranscriptComponent {
		}
	public interface Intron extends TranscriptComponent {
		}
	public interface CDS extends TranscriptComponent {
		}
	public interface UTR extends TranscriptComponent {
		}
	public interface UTR5 extends UTR {
		}
	public interface UTR3 extends UTR {
		}
	}
