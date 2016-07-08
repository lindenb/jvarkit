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
function GenomeBrowser()
	{
	this.minFontSize = 7;
	this.maxChromLength=100;
	this.useClip = false;
	this.width=1000;
	this.featureHeight=20;
	this.chromStart=0;
	this.chromEnd=0;
	this.bamFile='rf.bam';
	this.refFile='rf.fa';
	this.minHDistance=5;
	this.spaceYbetweenFeatures=2;
	this.minArrowWidth=2;
	this.maxArrowWidth=5;
	this.interval= null;
	this.genomicSequence = null;
	this.printReferenseSequence = true;
	this.printReadBases = true;
	this.printDeletions = false;
	this.printInsertions = false;
	this.hershey=new Hershey();
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


GenomeBrowser.prototype.paint=function(params)
	{

	var iter,rows=[];

	this.interval = null;

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
			height: margin_top+ rows.length*(this.spaceYbetweenFeatures+this.featureHeight)+this.spaceYbetweenFeatures
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
			ctx.lineWidth=0.5;
			ctx.strokeStyle=(x%vticks==0?"black":"gray");
			ctx.beginPath();
			ctx.moveTo(this.baseToPixel(x), 0);
			ctx.lineTo(this.baseToPixel(x), imageSize.height);
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
	var y=margin_top+this.spaceYbetweenFeatures;
	/** loop over each row */
	for(var rowY in rows)
		{

		var row=rows[rowY];
		for(x in row)
			{
			var rec=row[x];
			var x0=this.left(rec);

			var x1=this.right(rec);

			var y0=y;
			var y1=y0+this.featureHeight;
			var midY=y0+this.featureHeight/2.0;

			/* draw horizontal line */
			ctx.strokeStyle="black";

			ctx.moveTo(x0,midY);
			ctx.lineTo(x1,midY);
			ctx.stroke();



			/* draw record shape */
			ctx.beginPath();



			if(x1-x0 < this.minArrowWidth)
				{
				ctx.moveTo(x0,y0);
				ctx.lineTo(x1,y0);
				ctx.lineTo(x1,y1);
				ctx.lineTo(x0,y1);
				ctx.lineTo(x0,y0);
				}
			else
				{

				var arrow=Math.max(this.minArrowWidth,Math.min(this.maxArrowWidth, x1-x0));
				if(!rec.getReadNegativeStrandFlag())
					{
					ctx.moveTo(x0, y0);
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
			if(!rec.isReadPairedFlag() || rec.isProperPairFlag())
				{
				var grd=ctx.createLinearGradient(x0,y0,x0,y1);
				grd.addColorStop(0,"gray");
				grd.addColorStop(0.5,"white");
				grd.addColorStop(1.0,"gray");

				ctx.fillStyle=grd;
				}
			else
				{
				ctx.fillStyle="white";
				}
			ctx.fill();


			ctx.lineWidth=1;
			ctx.strokeStyle=this.getReadColor(rec);
			ctx.stroke();

			var cigarinder=0;
			var cigar = rec.getCigar();
			if(cigar==null || cigar.isEmpty())
				{
				continue;
				}
			var refpos=rec.getAlignmentStart();
			var readpos=0;

			/* loop over all cigar */
			for(cigarindex=0; cigarindex< cigar.getNumElements(); ++cigarindex)
				{
				var ce=cigar.get(cigarindex);
				var k=0;
				var mutW= this.baseToPixel(refpos+1) - this.baseToPixel(refpos);

				/* loop over this cigar-element */
				for(k=0;k< ce.getLength() ;++k)
					{
					var next_readpos=readpos;
					var next_refpos=refpos;


					ctx.lineWidth=1;
					ctx.strokeStyle="pink";
					switch(ce.getOperator().name)
						{
						case 'P':break;
						case 'N':
						case 'D':
							{
							if( this.printDeletions ) {
								ctx.fillStyle="red";
								ctx.beginPath();
								ctx.fillRect(
									this.baseToPixel(refpos),y0,
									mutW,y1-y0
									);
								ctx.closePath();
								}
							next_refpos++;
							break;
							}
						case 'I':
							{
							if( this.printInsertions ) {
								ctx.lineWidth=4;
								ctx.strokeStyle="yellow";
								ctx.beginPath();
								ctx.moveTo(this.baseToPixel(refpos),y0);
								ctx.lineTo(this.baseToPixel(refpos),y1);
								ctx.stroke();
								}
							next_readpos++;
							break;
							}
						case 'S':
						case 'H':
							{
							ctx.lineWidth=3;
							ctx.strokeStyle="green";
							ctx.beginPath();
							ctx.moveTo(this.baseToPixel(refpos),y0);
							ctx.lineTo(this.baseToPixel(refpos),y1);
							ctx.stroke();
							if(ce.getOperator().name=='S') next_readpos++;
							break;
							}
						case 'M':
						case 'X':
						case '=':
							{
							if( rec.hasReadString() && this.printReadBases && mutW>=this.minFontSize ) {

								var c1 = rec.getBaseAt(readpos);
								var c2= this.genomicSequence.charAt1(refpos);

								if(ce.getOperator().name == 'X' || (c1!='N' && c1!='n'  && c1.toUpperCase()!=c2.toUpperCase()))
									{
									ctx.lineWidth=1;
									ctx.strokeStyle="red";
									}
								else
									{
									ctx.lineWidth=0.3;
									ctx.strokeStyle="gray";
									}

								ctx.beginPath();
								this.hershey.paint(
									ctx,c1,
									this.baseToPixel(refpos),
									y0+2,
									mutW,
									(y1-y0)-4
									);
								ctx.stroke();
								}
							next_readpos++;
							next_refpos++;
							break;
							}
						default: console.log("unknown operator "+ce.getOperator().name);break;
						}

					readpos=next_readpos;
					refpos=next_refpos;

					}
				}

			}
		y+=this.featureHeight+this.spaceYbetweenFeatures;
		}

	/* draw frame */
	ctx.strokeStyle="black";
	ctx.rect(0,0,imageSize.width,imageSize.height);
	ctx.stroke();
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
