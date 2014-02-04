package com.github.lindenb.jvarkit.util.picard;

import java.util.HashSet;
import java.util.Set;

import net.sf.samtools.SAMRecord;

import com.github.lindenb.jvarkit.lang.Predicate;

public enum SamFlag
		{
		READ_PAIRED(0x1)
			{
			@Override
			public Predicate<SAMRecord> getAcceptFilter()
				{
				return new Predicate<SAMRecord>()
					{
					@Override
					public boolean apply(SAMRecord rec)
						{
						return rec.getReadPairedFlag();
						}
					};
				}	
			},
		READ_MAPPED_IN_PROPER_PAIR(0x2)
			{
			@Override
			public Predicate<SAMRecord> getAcceptFilter()
				{
				return new Predicate<SAMRecord>()
					{
					@Override
					public boolean apply(SAMRecord rec)
						{
						return rec.getReadPairedFlag() && rec.getProperPairFlag();
						}
					};
				}	
			},
		READ_UNMAPPED(0x4)
			{
			@Override
			public Predicate<SAMRecord> getAcceptFilter()
				{
				return new Predicate<SAMRecord>()
					{
					@Override
					public boolean apply(SAMRecord rec)
						{
						return rec.getReadUnmappedFlag();
						}
					};
				}	
			},
		MATE_UNMAPPED(0x8)
			{
			@Override
			public Predicate<SAMRecord> getAcceptFilter()
				{
				return new Predicate<SAMRecord>()
					{
					@Override
					public boolean apply(SAMRecord rec)
						{
						return rec.getReadPairedFlag() && rec.getMateUnmappedFlag();
						}
					};
				}	
			},
		READ_REVERSE_STRAND(0x10)
			{
			@Override
			public Predicate<SAMRecord> getAcceptFilter()
				{
				return new Predicate<SAMRecord>()
					{
					@Override
					public boolean apply(SAMRecord rec)
						{
						return rec.getReadNegativeStrandFlag();
						}
					};
				}	
			},
		MATE_REVERSE_STRAND(0x20)
			{
			@Override
			public Predicate<SAMRecord> getAcceptFilter()
				{
				return new Predicate<SAMRecord>()
					{
					@Override
					public boolean apply(SAMRecord rec)
						{
						return rec.getReadPairedFlag() && rec.getMateNegativeStrandFlag();
						}
					};
				}	
			},
		FIRST_IN_PAIR(0x40)
			{
			@Override
			public Predicate<SAMRecord> getAcceptFilter()
				{
				return new Predicate<SAMRecord>()
					{
					@Override
					public boolean apply(SAMRecord rec)
						{
						return rec.getReadPairedFlag() && rec.getFirstOfPairFlag();
						}
					};
				}	
			},
		SECOND_IN_PAIR(0x80)
			{
			@Override
			public Predicate<SAMRecord> getAcceptFilter()
				{
				return new Predicate<SAMRecord>()
					{
					@Override
					public boolean apply(SAMRecord rec)
						{
						return rec.getReadPairedFlag() && rec.getSecondOfPairFlag();
						}
					};
				}	
			},
		/*	the alignment is not primary (a read having split hits may have multiple primary alignment records).*/
		NOT_PRIMARY_ALIGNMENT(0x100)
			{
			@Override
			public Predicate<SAMRecord> getAcceptFilter()
				{
				return new Predicate<SAMRecord>()
					{
					@Override
					public boolean apply(SAMRecord rec)
						{
						return rec.getNotPrimaryAlignmentFlag();
						}
					};
				}	
			},
		READ_FAILS_VENDOR_QUALITY_CHECK(0x200)
			{
			@Override
			public Predicate<SAMRecord> getAcceptFilter()
				{
				return new Predicate<SAMRecord>()
					{
					@Override
					public boolean apply(SAMRecord rec)
						{
						return rec.getReadFailsVendorQualityCheckFlag();
						}
					};
				}	
			},
		READ_IS_DUPLICATE(0x400)
			{
			@Override
			public Predicate<SAMRecord> getAcceptFilter()
				{
				return new Predicate<SAMRecord>()
					{
					@Override
					public boolean apply(SAMRecord rec)
						{
						return rec.getDuplicateReadFlag();
						}
					};
				}	
			},
		SUPPLEMENTARY_ALIGNMENT(0x800)
			{
			@Override
			public Predicate<SAMRecord> getAcceptFilter()
				{
				return new Predicate<SAMRecord>()
					{
					@Override
					public boolean apply(SAMRecord rec)
						{
						return rec.getSupplementaryAlignmentFlag();
						}
					};
				}	
			}
		;
	private int flag;
	SamFlag(int flag)
		{
		this.flag=flag;
		}
	
	public int getFlag()
		{
		return flag;
		}

	public String getLabel()
		{
		return name().toLowerCase().replace('_', ' ');
		}
	
	public static SamFlag valueOf(int flg)
		{
		for(SamFlag f:values())
			{
			if(flg==f.getFlag()) return f;
			}
		return null;
		}
	
	public boolean isSet(int flag)
		{
		return (getFlag() & flag)!=0; 
		}
	
	public boolean isUnset(int flag)
		{
		return !isSet(flag);
		}
	
	public static Set<SamFlag> getFlags(int flag)
		{
		Set<SamFlag> set=new HashSet<SamFlag>();
		for(SamFlag f:values())
			{
			if(f.isSet(flag)) set.add(f);
			}
		return set;
		}
	
	public abstract Predicate<SAMRecord> getAcceptFilter();
	}
