/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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


/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
function Interval()
	{
	if(arguments.length==3)
		{
		this.contig = arguments[0];
		this.start = arguments[1];
		this.end = arguments[2];
		}
	else
		{
		var s= arguments[0];
		var colon=s.indexOf(':');
		if(colon<1)
			{
			throw "bad colon "+s;
			}
		var hyphen=s.indexOf('-');
		if(hyphen<colon)
			{
			throw "bad hyphen "+s;
			}
		this.contig = s.substr(0,colon);
		this.start = parseInt( s.substring(colon+1,hyphen).replace(/[,]/g, '') );
		this.end= parseInt( s.substring(hyphen+1).replace(/[,]/g, '') );
		}
	}
Interval.prototype.getContig = function() { return this.contig;}
Interval.prototype.getChrom = function() { return this.getContig();}
Interval.prototype.getStart = function() { return this.start;}
Interval.prototype.getEnd = function() { return this.end;}
Interval.prototype.distance = function() { return 1+(this.getEnd()-this.getStart());}
Interval.prototype.overlap = function(p) {
		return this.getContig()==p.getContig() && !(this.getEnd()<p.getStart() || p.getEnd()<this.getStart());
		}
Interval.prototype.toString = function() { return this.getContig()+":"+this.getStart()+"-"+this.getEnd();}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
function SamSequenceRecord(name,length,idx) {
	this.mSequenceName = name;
	this.mSequenceLength = length;
	this.mSequenceIndex = idx;
	}
SamSequenceRecord.prototype.getSequenceName =  function() { return this.mSequenceName; }
SamSequenceRecord.prototype.getSequenceLength =  function() { return this.mSequenceLength; }
SamSequenceRecord.prototype.getSequenceIndex =  function() { return this.mSequenceIndex; }
SamSequenceRecord.prototype.toString = function() { return this.getSequenceName()+":"+this.getSequenceLength();}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
function ReferenceSequence(name,index1,bases) {
	this.name = name;
	this.contigIndex1 = index1;
	this.bases = bases;
	}
ReferenceSequence.prototype.getSequenceName =  function() { return this.getContig(); }
ReferenceSequence.prototype.getName =  function() { return this.getContig(); }
ReferenceSequence.prototype.getContig =  function() { return this.name; }
ReferenceSequence.prototype.charAt = function(idx1)
			{
			if(this.bases==null) return 'N';
			if(idx1<this.contigIndex1) return 'N';
			idx1 -= this.contigIndex1;
			if(idx1>= this.bases.length ) return 'N';
			return this.str[idx1];
			}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
function CigarOperator(cread,cref,name)
	{
	this.name=name;
	this.cref = cref;
	this.cread = cread;
	}

CigarOperator.prototype.toString = function() { return this.name;}
CigarOperator.prototype.consumesReadBases = function() { return this.cread;}
CigarOperator.prototype.consumesReferenceBases = function() { return this.cref;}
CigarOperator.M = new CigarOperator(true, true,   'M');
CigarOperator.I = new CigarOperator(true, false,  'I');
CigarOperator.D = new CigarOperator(false, true,  'D');
CigarOperator.N = new CigarOperator(false, true,  'N');
CigarOperator.S = new CigarOperator(true, false,  'S');
CigarOperator.H = new CigarOperator(false, false, 'H');
CigarOperator.P = new CigarOperator(false, false, 'P');
CigarOperator.EQ = new CigarOperator(true, true,  '=');
CigarOperator.X = new CigarOperator(true, true,   'X');
CigarOperator.values = [
	CigarOperator.M,CigarOperator.I,CigarOperator.D,
	CigarOperator.S,CigarOperator.H,CigarOperator.P,
	CigarOperator.EQ,CigarOperator.X
	];
CigarOperator.characterToEnum = function(b) {
	if(b == "=") return CigarOperator.EQ;
	for(var i in CigarOperator.values )	{
		if( CigarOperator.values[i].name == b) return CigarOperator.values[i];
		}
	throw "bad cigar operator "+b;
	}
CigarOperator.prototype.isClip = function() {
	return this.name=="S" || this.name=="H";
	}
CigarOperator.prototype.isIndel = function() {
	return this.name=="I" || this.name=="D";
	}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

function CigarElement(op,length)
	{
	this.op=op;
	this.length = length;
	}
CigarElement.prototype.toString = function() { return ""+this.length+this.op.name;}
CigarElement.prototype.getOperator = function() { return this.op;}
CigarElement.prototype.getLength = function() { return this.length;}
CigarElement.prototype.size = function() { return this.getLength();}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
function Cigar()
	{
	this.elements=[];
	if( arguments.length == 1 )
		{
		var str = arguments[0];
		var i=0;
		while(i < str.length )
			{
			var prev=i;
			while( i< str.length )
				{
				var c = str.substring(i,i+1);
				if( c < '0' || c > '9') break;
				++i;
				}
			var len = parseInt(str.substring(prev,i) );
			var op =  CigarOperator.characterToEnum( str.substring(i,i+1) );
			this.elements.push( new CigarElement(op,len) );
			++i;
			}
		}
	}
Cigar.prototype.size = function() { return this.getNumElements();}
Cigar.prototype.getNumElements = function() { return this.elements.length;}
Cigar.prototype.isEmpty = function() { return this.getNumElements() == 0 ;}
Cigar.prototype.get = function(idx) { return this.elements[idx];}
Cigar.prototype.toString = function()
	{
	var s="";
	for(var i=0;i < this.getNumElements();++i) {
		s+= this.get(i).toString();
		}
	return s;
	};


Cigar.prototype.getUnclippedStart = function(pos) {
		var i=0;
        for (i=0;i< this.getNumElements();++i) {
        	if(this.get(i).getOperator().isClip()) {
                pos -= this.get(i).getLength();
            } else {
                break;
            }
        }
        return pos;
    };
    
Cigar.prototype.getUnclippedEnd = function(pos) {
        var i=0;
        for (i = this.getNumElements() - 1; i >= 0; --i) {
            if(this.get(i).getOperator().isClip()) {
                pos += this.get(i).getLength();
            } else {
                break;
            }
        }
        return pos;
    };

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/


function SamRecord()
	{
	this.name=null;
	this.flags=0;
	this.contig = null;
	this.pos = SamRecord.NO_ALIGNMENT_START ;
	this.mapq = SamRecord.NO_MAPPING_QUALITY ;
	this.cigarStr = null;
	this.cigar = null;
	this.alignEnd = null;
	this.sequence = SamRecord.NULL_SEQUENCE_STRING;
	this.qualities = null;
	
	if(arguments.length == 1) 
		{
		var read = arguments[0];
		if("name" in read) this.name = read.name;
		if("flag" in read) this.flags = read.flag;
		if("mapq" in read) this.mapq = read.mapq;
		if("pos" in read) this.pos = read.pos;
		if("sequence" in read) this.sequence = read.sequence;
		if("qualities" in read) this.qualities = read.qualities;
		if("cigar" in read) this.cigarStr = read.cigar;
		}
	
	};
	
SamRecord.NO_MAPPING_QUALITY = 0;
SamRecord.NO_ALIGNMENT_START = 0;
SamRecord.NULL_SEQUENCE_STRING = "*";
SamRecord.NO_ALIGNMENT_CIGAR = "*";

SamRecord.prototype.getReadName=function()
	{
	return this.name==null?"":this.name;
	};

SamRecord.prototype.setReadName=function(s)
	{
	this.name=s;
	return this;
	};

SamRecord.prototype.getName=function()
	{
	return this.getReadName();
	};

SamRecord.prototype.toString=function()
	{
	return this.getName()+" "+this.getReferenceName()+":"+this.getAlignmentStart()+"-"+this.getAlignmentEnd();
	};

SamRecord.prototype.setReferenceName=function(c)
	{
	this.contig = c;
	return this;
	};


SamRecord.prototype.getReferenceName=function()
	{
	return this.contig;
	};

SamRecord.prototype.setFlags=function(f)
	{
	this.flags = f;
	return this;
	};

SamRecord.prototype.getFlags=function()
	{
	return this.flags;
	};

SamRecord.prototype.isFlagSet=function(f)
	{
	return this.getFlags() & f;
	};

SamRecord.prototype.getReadPairedFlag=function()
	{
	return this.isReadPairedFlag();
	};

SamRecord.prototype.isReadPairedFlag=function()
	{
	return this.isFlagSet(0x1);
	};

SamRecord.prototype.getProperPairFlag=function()
	{
	return this.isProperPairFlag();
	};

SamRecord.prototype.isProperPairFlag=function()
	{
	if(!this.getReadPairedFlag() ) return false;
	return this.isFlagSet(0x2);
	};

SamRecord.prototype.isReadUnmappedFlag=function()
	{
	return this.isFlagSet(0x4);
	};

SamRecord.prototype.getReadUnmappedFlag=function()
	{
	return this.isReadUnmappedFlag();
	};

SamRecord.prototype.getMateUnmappedFlag=function()
	{
	if(!this.getReadPairedFlag() ) return false;
	return   this.isFlagSet(0x8);
	};
	
SamRecord.prototype.getReadNegativeStrandFlag=function()
	{
	return this.isFlagSet(0x10);
	};

SamRecord.prototype.getNotPrimaryAlignmentFlag=function()
	{
	return this.isFlagSet(0x100);
	};


SamRecord.prototype.isSupplementaryAlignmentFlag=function()
	{
	return this.isFlagSet(0x200);
	};	

SamRecord.prototype.getSupplementaryAlignmentFlag=function()
	{
	return this.isSupplementaryAlignmentFlag();
	};	

SamRecord.prototype.getDuplicateReadFlag=function()
	{
	return this.isFlagSet(0x400);
	};
	

SamRecord.prototype.getReadFailsVendorQualityCheckFlag=function()
	{
	return this.isFlagSet(0x800);
	};

SamRecord.prototype.setAlignmentStart=function(v)
	{
	this.pos = v;
	this.alignEnd = null;
	return this;
	};

SamRecord.prototype.getAlignmentStart=function()
	{
	return this.pos;
	};

SamRecord.prototype.getStart=function()
	{
	return this.getAlignmentStart();
	};

	
SamRecord.prototype.getAlignmentEnd=function()
	{
	if(this.getReadUnmappedFlag()) return SamRecord.NO_ALIGNMENT_START;
	if(this.alignEnd==null)
		{
		var cigar = this.getCigar();
		if( cigar==null) throw "cigar == null";
		this.alignEnd=this.getAlignmentStart();
		for (var i =0;i< cigar.getNumElements();++i)
		  	{
		  	var ce = cigar.get(i);
		  	if( ce.getOperator().consumesReferenceBases() )
		  		{
		  		this.alignEnd += ce.getLength() ;
			 	}
            }
        }
    return this.alignEnd;
	};

SamRecord.prototype.getEnd = function()
	{
	return this.getAlignmentEnd();
	};
	
	
SamRecord.prototype.getUnclippedStart=function() {
        return this.getCigar().getUnclippedStart(this.getAlignmentStart());
    }

SamRecord.prototype.getUnclippedEnd=function() {
        return this.getCigar().getUnclippedEnd(this.getAlignmentEnd());
    };
	
SamRecord.prototype.getCigarString=function() {
	return this.cigarStr;
	}

SamRecord.prototype.getCigar=function()
	{
	if(this.cigar!=null) return this.cigar;
	if(this.isReadUnmappedFlag()) return null;
	var s = this.getCigarString();
	
	if(s  == null || s == SamRecord.NO_ALIGNMENT_CIGAR ) return null;
	
	this.cigar = new Cigar(this.getCigarString());
	return this.cigar;
	};

SamRecord.prototype.hasReadString = function() {
	var s= this.getReadString();
	return s!=null && s != SamRecord.NULL_SEQUENCE_STRING &&  s.length>0;
	};

SamRecord.prototype.getReadString = function() {
	return this.sequence;
	};

SamRecord.prototype.setReadString = function(v) {
	this.sequence = v;
	return this;
	};

SamRecord.prototype.getBaseAt  = function(idx)
	{
	if(!this.hasReadString()) return '?';
	var s = this.getReadString();
	if(idx<0 || idx>= s.length) return '?';
    return s.charAt(idx);
	};

/*
var x = new SamRecord({"name":"HWI-1KL149:110:C8HU4ACXX:8:1216:15251:82595","flag":83,"ref":"1","pos":1647920,"mapq":60,"cigar":"49S2M49S","len":-153,"materef":"1","matepos":1647865,"sequence":"AATGAGAAATAAAGTGTCATGCAAAGAAACCTCACTTCAAAAATTTCACATGAAGCCGGGCACGGAGGCTTATGCCTGTAATCCTAGCACTTTGGGAGGC","qualities":"CBBEGFGGGFGGIHGHFIGHIIGGFGFFEFGEHDGEDGDDDDEDDDFBFEEEDEECBCCDEBABDDBDECCDDDBEDDCCDCBECDDFAFEEEFECC?BB"});
print(x.getCigar());
print(x.getUnclippedStart());
print(x.getUnclippedEnd());
*/
