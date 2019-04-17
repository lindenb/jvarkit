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
package com.github.lindenb.jvarkit.util.bio.samfilter;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.function.Predicate;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.ParameterException;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.RuntimeIOException;
/**

Note to self:

java -cp ~/package/javacc/javacc.jar javacc -OUTPUT_DIRECTORY=src/main/java/com/github/lindenb/jvarkit/util/bio/samfilter -JDK_VERSION=1.8 src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 

*/
public class SamRecordFilterFactory  implements IStringConverter<SamRecordFilter> {

private static final Logger LOG = Logger.build(SamRecordFilterFactory.class).make();

public static final String FILTER_DESCRIPTION = "A filter expression. Reads matching the expression will be filtered-out. Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj for a complete syntax. ";
public static final String DEFAULT_FILTER = "mapqlt(1) || MapQUnavailable() || Duplicate() || FailsVendorQuality() || NotPrimaryAlignment() || SupplementaryAlignment()";
public static final String DEFAULT_OPT = "--samFilter";

/** a SamRecordFilter accepting any SAMRecord */
public static final SamRecordFilter  ACCEPT_ALL = new SamRecordFilter() {
    @Override
    public boolean filterOut(final SAMRecord first, SAMRecord second) {
            return false;
    }

    @Override
    public boolean filterOut(SAMRecord record) {
            return false;
    }
    @Override
    public String toString() {
            return "Accept All/ Filter out nothing";
    }
};

/** get Default filter using a static way */
public static SamRecordFilter getDefault() {
     return new SamRecordFilterFactory().buildDefault();
    }



public SamRecordFilter buildDefault() {
    try {
        return build(DEFAULT_FILTER);
        }
    catch(final ParseException err)
        {
        throw new IllegalArgumentException(err);
        }
    }

@Override
/** if query is empty/null , return ACCEPT_ALL */
public SamRecordFilter convert(final String query) {
     if(StringUtils.isBlank(query)) return ACCEPT_ALL;
     try {
        return build(query);
     	}
     catch (final ParseException e) {
         throw new ParameterException("Cannot parse samFilter : "+query,e);
         }
	}


Predicate<SAMRecord> overlapBed(final String fname) {
    final Path bedFile = Paths.get(fname);
    IOUtil.assertFileIsReadable(bedFile);
    final BedLineCodec codec = new BedLineCodec();
    final IntervalTreeMap<Boolean> intervals = new IntervalTreeMap<Boolean>();
  
    try (BufferedReader r = IOUtils.openPathForBufferedReading(bedFile)) {
    	String line;
    	while((line=r.readLine())!=null) {
    		final BedLine bedline = codec.decode(line);
    		if(bedline==null)  continue;
    		intervals.put(bedline.toInterval(),Boolean.TRUE);
    		}
        return new Predicate<SAMRecord>() {
			@Override
			public boolean test(final SAMRecord t) {
				return !t.getReadUnmappedFlag() &&  intervals.containsOverlapping(t)
						;
			}
		};
    } catch(final IOException err) {
    	LOG.error(err);
    	throw new RuntimeIOException(err);
    }
  
}


Predicate<SAMRecord> duplicateFilter() { return rec->rec.getDuplicateReadFlag();}

Predicate<SAMRecord> unmappedFilter()  { return rec->rec.getReadUnmappedFlag();}

Predicate<SAMRecord> mappedFilter() { return rec->!rec.getReadUnmappedFlag();}

Predicate<SAMRecord> failsVendorQuality() { return rec-> rec.getReadFailsVendorQualityCheckFlag();}

Predicate<SAMRecord> readPaired() {
	return new Predicate<SAMRecord>() {
	 	@Override public boolean test(final SAMRecord rec) { return rec.getReadPairedFlag();}
	}; }
Predicate<SAMRecord> mateUnmapped() {
	return new Predicate<SAMRecord>() {
	 	@Override public boolean test(final SAMRecord rec) { return rec.getMateUnmappedFlag();}
	}; }
Predicate<SAMRecord> samFlag(final int flg) {
	return new Predicate<SAMRecord>() {
	 	@Override public boolean test(final SAMRecord rec) { return (rec.getFlags() & flg) != 0;}
	}; }
Predicate<SAMRecord> sample(final String s) {
    return new Predicate<SAMRecord>() {
            @Override public boolean test(final SAMRecord rec) { SAMReadGroupRecord rg=rec.getReadGroup(); return rg!=null && s.equals(rg.getSample());}
    }; }

Predicate<SAMRecord> group(final String s) {
    return new Predicate<SAMRecord>() {
            @Override public boolean test(final SAMRecord rec) { SAMReadGroupRecord rg=rec.getReadGroup(); return rg!=null && s.equals(rg.getId());}
    }; }
Predicate<SAMRecord> notPrimaryAlignmentFlag() {
return new Predicate<SAMRecord>() {
 	@Override public boolean test(final SAMRecord rec) { return rec.isSecondaryAlignment();}
}; }		

Predicate<SAMRecord> supplementaryAlignmentFlag() {
return new Predicate<SAMRecord>() {
 	@Override public boolean test(final SAMRecord rec) { return rec.getSupplementaryAlignmentFlag();}
}; }		

Predicate<SAMRecord> readClipped() {
return new Predicate<SAMRecord>() {
        @Override public boolean test(final SAMRecord rec) {
        if(rec.getReadUnmappedFlag()) return false;
        final Cigar c= rec.getCigar();
        if(c==null || c.isEmpty()) return false;
        return c.isClipped();
        	
        }
}; }

Predicate<SAMRecord>  mapqUnavailable() {
return new Predicate<SAMRecord>() {
	@Override public boolean test(final SAMRecord rec) { return  (rec.getMappingQuality() == SAMRecord.NO_MAPPING_QUALITY);}
}; }		



Predicate<SAMRecord>  hasFlag(final int flg) {
return new Predicate<SAMRecord>() {
	@Override public boolean test(final SAMRecord rec) { return   (rec.getFlags() & flg) != 0;}
}; }	


Predicate<SAMRecord>  discordant() {
    return rec->
    	  		rec.getReadPairedFlag() &&
    			!rec.getReadUnmappedFlag() &&
    			!rec.getMateUnmappedFlag() &&
    			rec.getReferenceIndex()!=rec.getMateReferenceIndex();
    }


Predicate<SAMRecord>  mapqLowerThan(final int mapq) { return rec-> rec.getMappingQuality() < mapq ;}


/** parse predicate returning **true** if the record should be **REJECTED**/
public Predicate<SAMRecord> parseRejectPredicate(final String query) throws ParseException {
	 if(StringUtils.isBlank(query)) return SR->false;// never reject anything
	 
	 try(Reader r= new StringReader(query)) {
		final SamFilterParser sfp =new SamFilterParser(r);
		return  sfp.anyNode();
	 	}
     catch (final Exception e) {
         throw new ParameterException("Cannot parse samFilter : "+query,e);
         }
	}

private SamRecordFilter build(final String expr) throws ParseException {
	    final Predicate<SAMRecord> pred = parseRejectPredicate(expr);
	    
	    return new SamRecordFilter() {
            @Override
            public boolean filterOut(final SAMRecord first, final SAMRecord second) {
                    throw new IllegalStateException("SamRecordFilter.filterOut(a,b); <- shouldn't happen");
            }

            @Override
            public boolean filterOut(final SAMRecord record) {
                    return pred.test(record);
            }
            @Override
            public String toString() {
                return expr;
                }
	    	};
		}


}
