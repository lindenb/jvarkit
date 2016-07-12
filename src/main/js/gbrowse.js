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

function Consensus()
	{
	this.depth=0;
	this.total=0;
	this.bases={};
	}
Consensus.prototype.watch = function(base) {
	if(!(base in this.bases)) {
		this.bases[base]=1;
		}
	else
		{
		this.bases[base]++;
		}
	this.total++;
	if(base!='-') this.depth++;
	};
	


function GenomeBrowser()
	{
	this.minFontSize = 7;
	this.maxChromLength=10000;
	this.useClip = false;
	this.width=1000;
	this.featureHeight=20;
	this.chromStart=0;
	this.chromEnd=0;
	this.minHDistance=5;
	this.spaceYbetweenFeatures=2;
	this.minArrowWidth=2;
	this.maxArrowWidth=5;
	this.interval= null;
	this.genomicSequence = null;
	this.coverageHeight=100;
	this.consensusHeight=this.featureHeight,
	this.printReferenseSequence = true;
	this.printReadBases = true;
	this.printReadName = false;
	this.printDeletions = false;
	this.printInsertions = false;
	this.expandInsertions = false;
	this.expandDeletions = true;
	this.hershey=new Hershey();
	this.base2consensus={};
	}

/* pretty print a number, with commas */
GenomeBrowser.prototype.niceNumber=function(n)
	{
	var y = "";
	var arr = parseInt(n).toString().split("");
	for(var i=0; i<arr.length; i++)
		{
	    	y += arr[i];
	    	if((arr.length-i-1)%3==0 && i<arr.length-1) y += ",";
		}
	return y;
	}

/* convert position to pixel */
GenomeBrowser.prototype.baseToPixel=function(pos)
	{
	return ((pos-this.interval.start)/(this.interval.end-this.interval.start))*this.width;
	};


GenomeBrowser.prototype.readStart=function(rec)
	{
	if(this.useClip)
		{
		return rec.getUnclippedStart();
		}
	else
		{
		return rec.getAlignmentStart();
		}
	}

GenomeBrowser.prototype.readEnd=function(rec)
	{
	if(this.useClip)
		{
		return rec.getUnclippedEnd();
		}
	else
		{
		return rec.getAlignmentEnd();
		}
	}

GenomeBrowser.prototype.left=function(rec)
	{
	return this.baseToPixel(this.readStart(rec));
	};

GenomeBrowser.prototype.right=function(rec)
	{
	return this.baseToPixel(this.readEnd(rec));
	};

GenomeBrowser.prototype.getReadColor=function(rec)
	{
	if(rec.getMateUnmappedFlag()) return "green";
	if(rec.getNotPrimaryAlignmentFlag()) return "pink";
	if(rec.getDuplicateReadFlag()) return "fuchsia";
	if(rec.getSupplementaryAlignmentFlag()) return "skyblue";
	if(rec.getReadFailsVendorQualityCheckFlag()) return "khaki";
	return "black";
	}

GenomeBrowser.prototype.paintCigarComponent = function(ctx,r) {
	if( r.next_refpos < this.interval.start) return;
	if( r.refpos > this.interval.end) return;
	
	
	var mutW= this.baseToPixel(r.refpos+1) - this.baseToPixel(r.refpos);
	var x0 = this.baseToPixel(r.refpos);
	var x1 = this.baseToPixel(r.next_refpos);
	var y0 = r.y;
	var y1 = y0 + this.featureHeight-2;
	var midY= y0 + (this.featureHeight-2)/2.0;
	
	if(r.ce.getOperator().isOneOf("DN"))
		{
		y0 = midY;
		y1 = midY;
		}
	
	/* draw horizontal line */
	ctx.globalAlpha = 1;
	ctx.lineWidth=1;
	ctx.strokeStyle="black";
	ctx.fillStyle="black";
	ctx.beginPath();
	ctx.moveTo(x0,midY);
	ctx.lineTo(x1,midY);
	ctx.closePath();
	ctx.stroke();

	
	
	/* draw record shape */
	ctx.beginPath();

	if(r.arrow==0 || mutW < this.minArrowWidth)
		{
		ctx.moveTo(x0,y0);
		ctx.lineTo(x1,y0);
		ctx.lineTo(x1,y1);
		ctx.lineTo(x0,y1);
		ctx.lineTo(x0,y0);
		}
	else
		{
		var arrow=Math.max(this.minArrowWidth,Math.min(this.maxArrowWidth, x1 - x0));
		if(r.arrow==1)
			{
			ctx.moveTo(x0,y0);
			ctx.lineTo(x1-arrow,y0);
			ctx.lineTo(x1,midY);
			ctx.lineTo(x1-arrow,y1);
			ctx.lineTo(x0,y1);
			}
		else
			{
			ctx.moveTo(x0+arrow, y0);
			ctx.lineTo(x0,midY);
			ctx.lineTo(x0+arrow,y1);
			ctx.lineTo(x1,y1);
			ctx.lineTo(x1,y0);
			}
		}
	ctx.closePath();
	
	
	
	if(r.ce.getOperator().isOneOf("M=") &&
			r.rec.isReadPairedFlag() &&
			r.rec.isProperPairFlag() && 
			!r.rec.hasDiscordantContigs() &&
			!r.rec.isNotPrimaryAlignmentFlag() &&
			!r.rec.isSupplementaryAlignmentFlag()
			)
		{
		var grd=ctx.createLinearGradient(x0,y0,x0,y1);
		grd.addColorStop(0,"gray");
		grd.addColorStop(0.5,"white");
		grd.addColorStop(1.0,"gray");
		ctx.fillStyle=grd;
		
		}
	else if(r.ce.getOperator().name=="S")
		{
		ctx.fillStyle="pink";
		}
	else if(r.ce.getOperator().name=="I")
		{
		ctx.fillStyle="orange";
		}
	else if(r.rec.hasDiscordantContigs())
		{
		ctx.fillStyle="blue";
		}
	else
		{
		ctx.fillStyle="white";
		}
	
	
	ctx.lineWidth=1;
	ctx.strokeStyle="black";
	if(r.ce.getOperator().isOneOf("IDN"))
		{
		ctx.strokeStyle="red";
		ctx.lineWidth=5;
		}
	
	ctx.fill();
	ctx.stroke();
	
	/* update consensus for deletions */
	if(r.ce.getOperator().isOneOf("DN"))
		{
		for(var x=0 ; x < r.ce.getLength();++x) {
			if( r.refpos + x < this.interval.start) continue;
			if( r.refpos + x > this.interval.end) break;
			this.consensus(r.refpos+x,'-');
			}
		}
	
	
	var printBases = this.printReadBases && mutW>=this.minFontSize && r.readpos < r.next_readpos &&
					 ( r.ce.getOperator().isOneOf("X=MS") || (r.ce.getOperator().name == 'I' && this.expandInsertions));
			
	ctx.lineWidth=1;
	ctx.save();
	ctx.clip();
	
		
		for(var x=0 ; x < r.ce.getLength();++x) {
			if( r.refpos + x < this.interval.start) continue;
			if( r.refpos + x > this.interval.end) break;
		
			var c1 =r.rec.getBaseAt(r.readpos + x );
			var c2=  this.genomicSequence.charAt1(r.refpos + x);
			if(r.ce.getOperator().name == 'X' || (c1!='N' && c1!='n' && c2!='N' && c2!='n'  && c1.toUpperCase()!=c2.toUpperCase()))
				{
				ctx.beginPath();
				ctx.fillStyle="coral";
				ctx.strokeStyle="orange";
				ctx.rect(this.baseToPixel(r.refpos+ x),y0,mutW,(y1-y0));
				ctx.fill();
				ctx.stroke();
				}

			this.consensus(r.refpos+x,c1);
			
			if( !printBases ) continue;
			ctx.strokeStyle=this.base2color(c1);
			
			if( this.printReadName )
				{
				var readName = r.rec.getReadName();
				c1 = r.readpos +x<readName.length? readName.charAt(r.readpos +x):' ';
				ctx.strokeStyle="black";
				}
		
			
			
			ctx.beginPath();
			this.hershey.paint(
				ctx,c1,
				this.baseToPixel(r.refpos+ x),
				y0+2,
				mutW,
				(y1-y0)-4
				);
			//ctx.closePath(); NO, would could 'C' -> 'O'
			ctx.stroke();
			}
		
	ctx.restore();
	}

GenomeBrowser.prototype.consensus = function(pos,base) {
	var posstr=""+pos,c;
	if(!(posstr in this.base2consensus)) {
		c = new Consensus();
		this.base2consensus[posstr] = c;
		}
	else
		{
		c= this.base2consensus[posstr] ;
		}
	c.watch(base);
	}

GenomeBrowser.prototype.paint=function(params)
	{

	var iter,rows=[];

	this.interval = null;
	this.base2consensus={};
	
	if( !("canvasid" in params) )
		{
		throw "no @canvasid in params";
		}

	if( !("reads" in params) )
		{
		throw "no @reads in params";
		}


	if( !("reference" in params) )
		{
		this.genomicSequence = new ReferenceSequence("Un",1,"");
		}
	else
		{
		this.genomicSequence = params.reference ;
		}

	if( "interval" in params)
		{
		this.interval = params.interval;
		}
	else /* no interval provided, use reads to gather interval */
		{
		var chr=null;
		var B = null;
		var E = null;
		for(iter in params.reads )
			{
			var rec= params.reads[iter];
			if(rec.isReadUnmappedFlag()) continue;
			if( chr==null ) chr= rec.getReferenceName();
			else if( chr !=  rec.getReferenceName()) continue;
			var pos  =  this.readStart(rec);
			B = (B==null ? pos : pos < B ? pos : B );
			pos  =  this.readEnd(rec);
			E = (E==null ? pos : pos > E ? pos : E );
			}
		this.interval = new Interval(
			(chr==null?"Undefined":chr),
			(B==null?0:B),
			(E==null?0:E)
			);
		}


	if(this.interval.start > this.interval.end )
		{
		throw "bad interval "+this.interval;
		}
	if(this.interval.start + this.maxChromLength < this.interval.end )
		{
		this.interval.end = this.interval.start +  this.maxChromLength;
		}



	/** pileup reads */
	for(iter in params.reads )
		{

		var rec= params.reads[iter];

		if(rec.isReadUnmappedFlag()) continue;

		if( this.readEnd(rec) < this.interval.getStart() ) continue;
		if( this.readStart(rec) > this.interval.getEnd() ) continue;

		var y=0;
		for(y=0;y< rows.length;++y)
			{
			var row=rows[y];
			var last=row[row.length-1];
			if(this.right(last)+ this.minHDistance > this.left(rec)) continue;

			row.push(rec);
			rec=null;
			break;
			}
		if(rec!=null)
			{
			var row=[];
			row.push(rec);
			rows.push(row);
			}
		}




	var ruler_height = this.niceNumber( this.interval.getEnd() ).length*20;

	var refw= this.width/(1.0*this.interval.distance() ); if(refw<1.0) refw=1;refw=parseInt(refw);

	var margin_top=10+(refw*2)+ruler_height;
	var imageSize=	{
			width:  this.width,
			height: margin_top+ this.coverageHeight + this.consensusHeight + rows.length*(this.spaceYbetweenFeatures+this.featureHeight)+this.spaceYbetweenFeatures
			};

	var canvasdoc=document.getElementById(params.canvasid);
	if( canvasdoc == null) throw "Cannot get canvas id "+ params.canvasid;
	canvasdoc.setAttribute("width",imageSize.width);
	canvasdoc.setAttribute("height",imageSize.height);
	if(!canvasdoc.getContext)return;
	var ctx = canvasdoc.getContext('2d');
	if(ctx==null) return;



	/* draw reference sequence */
	var vticks = parseInt(Math.ceil((this.interval.distance()/10.0)));
	if( vticks< 1)   vticks =1;
	vticks = Math.pow(10,Math.floor(Math.log10(vticks)));



	for(var x= this.interval.start;
		x<=  this.interval.end;
		++x)
		{

		var oneBaseWidth = this.baseToPixel(x+1)-this.baseToPixel(x);


		if( oneBaseWidth > 4.0 || x%vticks==0)
			{
			//draw vertical line
			ctx.beginPath();
			ctx.lineWidth=(x%vticks==0?0.5:0.2);
			ctx.strokeStyle=(x%vticks==0?"black":"gray");
			ctx.moveTo( this.baseToPixel(x), 0);
			ctx.lineTo( this.baseToPixel(x), imageSize.height);
			ctx.closePath();
			ctx.stroke();
			}



		//draw base position
		if((x)%vticks==0 || oneBaseWidth >= this.minFontSize)
			{
			ctx.strokeStyle="black";
			var xStr=this.niceNumber(x);
			ctx.save();

			var numHeight = ((x)%vticks==0 && oneBaseWidth < this.minFontSize?7:oneBaseWidth-2);

			ctx.translate(this.baseToPixel(x)+numHeight,0);
			ctx.rotate(Math.PI/2.0);

			ctx.beginPath();
			this.hershey.paint(ctx,
					xStr,
					0,
					0,
					ruler_height,
					numHeight
					);
			ctx.stroke();
			ctx.restore();
			}

		if( this.printReferenseSequence && oneBaseWidth >= this.minFontSize)
			{
			//paint genomic sequence
			var c=this.genomicSequence.charAt1(x);
			ctx.lineWidth = 1;
			ctx.strokeStyle= (this.base2color(c)) ;
			ctx.beginPath();
			this.hershey.paint(ctx,
					c,
					this.baseToPixel(x)+1,
					ruler_height,
					oneBaseWidth-2,
					oneBaseWidth-2
					);
			ctx.stroke();
			}


		}

	ctx.lineWidth=1;
	var y=margin_top+ this.spaceYbetweenFeatures + this.coverageHeight + this.consensusHeight;
	/** loop over each row */
	for(var rowY in rows)
		{

		var row=rows[rowY];
		for(x in row)
			{
			var rec=row[x];
			var cigarindex=0;
			var cigar = rec.getCigar();
			if(cigar==null || cigar.isEmpty())
				{
				continue;
				}
			var refpos = rec.getUnclippedStart();
			var readpos = 0;
			var cigarcomponents=[];
			
			/* loop over all cigar */
			for(cigarindex=0; cigarindex< cigar.getNumElements(); ++cigarindex)
				{
				var ce=cigar.get(cigarindex);
				if( ce.getOperator().name == 'P') continue;
				var next_refpos=refpos;
				var next_readpos=readpos;
				
				var cigarcomponent={
					"y":y,
					"ce":ce,
					"rec":rec,
					"refpos":refpos,
					"readpos":readpos,
					"next_refpos": refpos,
					"next_readpos": readpos,
					"seq":"",
					"arrow":0
					};
				
				
				
				switch(ce.getOperator().name)
					{
					case 'P': break;
					case 'N':
					case 'D':
						{
						if(this.expandDeletions) {
							next_refpos += ce.getLength();
							}
						break;
						}
					case 'I':
						{
						if( this.expandInsertions )
							{
							next_refpos += ce.getLength();
							}
						next_readpos += ce.getLength();
						break;
						}
					case 'S':
						{
						next_refpos += ce.getLength();
						next_readpos += ce.getLength();
						break;
						}
					case 'H':
						{
						next_refpos += ce.getLength();
						break;
						}
					case 'M':
					case 'X':
					case '=':
						{
						next_refpos += ce.getLength();
						next_readpos += ce.getLength();
						break;
						}
					default: console.log("unknown operator "+ce.getOperator().name);break;
					}
				cigarcomponent.next_readpos = next_readpos;
				cigarcomponent.next_refpos = next_refpos;
				readpos = next_readpos;
				refpos = next_refpos;
				if( !this.useClip && ce.getOperator().isOneOf("SH")) continue;
				cigarcomponents.push( 	cigarcomponent );
				}
			if( rec.isReadNegativeStrandFlag() )
				{
				cigarcomponents[ 0 ].arrow = -1;
				}
			else
				{
				cigarcomponents[ cigarcomponents.length-1 ].arrow = 1;
				}
			//sort so, deletion are sorted first, insertions are last
			cigarcomponents.sort( function(r1, r2) {
				var i1=(r1.ce.hasOperatorIn("DN")?0:r1.ce.hasOperatorIn("I")?2:1);
				var i2=(r2.ce.hasOperatorIn("DN")?0:r2.ce.hasOperatorIn("I")?2:1);
				return i1-i2;
				});
			
			for(cigarindex in cigarcomponents)
				{
				this.paintCigarComponent(ctx,cigarcomponents[cigarindex]);
				}
			}
		y+=this.featureHeight+this.spaceYbetweenFeatures;
		}
	
	y= margin_top;
	for(var x= this.interval.start;
		x<=  this.interval.end;
		++x)
		{
		var base,c,w,x0 = this.baseToPixel(x);
		var oneBaseWidth = this.baseToPixel(x+1)-x0;
		if( oneBaseWidth < this.minFontSize ) break;
		var xstr = ""+x;
		if(!(xstr in this.base2consensus)) continue;
		var c= this.base2consensus[xstr];
		x0++;
		//loop over bases
		for(base in c.bases)
			{
			w = (c.bases[base]/(1.0*c.total))*(oneBaseWidth-2);
			ctx.strokeStyle= this.base2color(base) ;
			ctx.beginPath();
			this.hershey.paint(ctx,
					base,
					x0,
					y,
					w, oneBaseWidth-2
					);
			ctx.stroke();
			x0+=w;
			}
	
		}
		
		
	/* draw coverage */
	
	var maxdepth = 0.0;
	for(var x in this.base2consensus) {
		maxdepth= Math.max( maxdepth, this.base2consensus[x].depth );
		}
	
	y= margin_top+ this.spaceYbetweenFeatures + this.coverageHeight + this.consensusHeight;
	
	ctx.fillStyle="gray";
	ctx.strokeStyle="lightgray";
	for(var x= this.interval.start;
		maxdepth>0 && x<=  this.interval.end;
		++x)
		{
		var xstr= ""+x;
		var dp = ((( xstr in this.base2consensus ?this.base2consensus[xstr].depth : 0.0)*1.0 )/ maxdepth)*this.coverageHeight;
		
		var x0 = this.baseToPixel(x);
		var x1 = this.baseToPixel(x+1);
		if( (x1-x0) < this.minFontSize ) break;
		ctx.beginPath();
		ctx.moveTo(0,y);
		ctx.lineTo(  x0,y);
		ctx.lineTo(  x0, y-dp);
		ctx.lineTo(  x1, y-dp);
		ctx.lineTo(  x1, y);
		ctx.closePath();
		ctx.fill();
		ctx.stroke();
		}
	
	
	
	ctx.fillStyle="black";
	ctx.font = "12px serif";
	ctx.fillText("Max Depth:"+maxdepth, 10, y-12);
		

	/* draw frame */
	ctx.beginPath();
	ctx.strokeStyle="black";
	ctx.rect(0,0,imageSize.width,imageSize.height);
	ctx.stroke();
	
	
	
	this.base2consensus={};
	}



GenomeBrowser.prototype.base2color=function(c)
	{
	switch(c)
		{
		case 'n':case 'N': return "black";
		case 'a':case 'A': return "red";
		case 't':case 'T': return "green";
		case 'g':case 'G': return "yellow";
		case 'c':case 'C': return "blue";
		default: return "orange";
		}
	};
