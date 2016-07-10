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
function Hershey()
	{
	this.scalex=10.0;
	this.scaley=10.0;
	
	this.LETTERS=[
		"  9MWRMNV RRMVV RPSTS",
		" 16MWOMOV ROMSMUNUPSQ ROQSQURUUSVOV",
		" 11MXVNTMRMPNOPOSPURVTVVU",
		" 12MWOMOV ROMRMTNUPUSTURVOV",
		" 12MWOMOV ROMUM ROQSQ ROVUV",
		"  9MVOMOV ROMUM ROQSQ",
		" 15MXVNTMRMPNOPOSPURVTVVUVR RSRVR",
		"  9MWOMOV RUMUV ROQUQ",
		"  3PTRMRV",
		"  7NUSMSTRVPVOTOS",
		"  9MWOMOV RUMOS RQQUV",
		"  6MVOMOV ROVUV",
		" 12LXNMNV RNMRV RVMRV RVMVV",
		"  9MWOMOV ROMUV RUMUV",
		" 14MXRMPNOPOSPURVSVUUVSVPUNSMRM",
		" 10MWOMOV ROMSMUNUQSROR",
		" 17MXRMPNOPOSPURVSVUUVSVPUNSMRM RSTVW",
		" 13MWOMOV ROMSMUNUQSROR RRRUV",
		" 13MWUNSMQMONOOPPTRUSUUSVQVOU",
		"  6MWRMRV RNMVM",
		"  9MXOMOSPURVSVUUVSVM",
		"  6MWNMRV RVMRV",
		" 12LXNMPV RRMPV RRMTV RVMTV",
		"  6MWOMUV RUMOV",
		"  7MWNMRQRV RVMRQ",
		"  9MWUMOV ROMUM ROVUV"
		];
	this.DIGITS=[
		  " 12MWRMPNOPOSPURVTUUSUPTNRM",
		  "  4MWPORMRV",
		  "  9MWONQMSMUNUPTROVUV",
		  " 15MWONQMSMUNUPSQ RRQSQURUUSVQVOU",
		  "  7MWSMSV RSMNSVS",
		  " 14MWPMOQQPRPTQUSTURVQVOU RPMTM",
		  " 14MWTMRMPNOPOSPURVTUUSTQRPPQOS",
		  "  6MWUMQV ROMUM",
		  " 19MWQMONOPQQSQUPUNSMQM RQQOROUQVSVUUURSQ",
		  " 14MWUPTRRSPROPPNRMTNUPUSTURVPV"
		];
	 }

Hershey.MOVETO=0;
Hershey.LINETO=1;
	
	
	
Hershey.prototype.charToHersheyString=function(c)
		{
		var codeA = "A".charCodeAt(0);
		var code0 = "0".charCodeAt(0);
		if(	c.toUpperCase().charCodeAt(0)>= codeA && 
			c.toUpperCase().charCodeAt(0)<="Z".charCodeAt(0)
			)
			{
		
			var idx=c.toUpperCase().charCodeAt(0)-codeA;
			return this.LETTERS[idx];
			}
		if(	c.charCodeAt(0)>=code0 && 
			c.charCodeAt(0)<="9".charCodeAt(0)
			)
			{
			return this.DIGITS[c.charCodeAt(0) - code0];
			}
		switch(c)
			{
			case '.': return "  6PURURVSVSURU";//210
			case ',': return "  7PUSVRVRUSUSWRY";//211
			case ':': return " 12PURPRQSQSPRP RRURVSVSURU";//212
			case ';': return " 13PURPRQSQSPRP RSVRVRUSUSWRY";//213
			case '!': return " 12PURMRR RSMSR RRURVSVSURU";//214
			case '?': return " 17NWPNRMSMUNUPRQRRSRSQUP RRURVSVSURU";//215
			case '\'':return "  3PTRMRQ";//216
			case '\"':return "  6NVPMPQ RTMTQ";//217
			case '/': return "  3MWVLNW";//220
			case '(': return "  7OVTLRNQPQSRUTW";//221
			case ')': return "  7NUPLRNSPSSRUPW";//222
			case '|': return "  3PTRLRW";//223
			case '#': return " 12MXRLPW RULSW ROPVP ROSVS";//233
			case '*': return "  9JZRLRX RMOWU RWOMU";//728
			case '=': return "  6LXNPVP RNTVT";//226
			case '-': return "  3KYKRYR";//806
			case '_': return "  3JZJZZZ";//998
			case '[': return " 12MWPHP\\ RQHQ\\ RPHUH RP\\U\\";//1223
			case ']': return " 12MWSHS\\ RTHT\\ ROHTH RO\\T\\";//1224
			case '{': return " 38LWSHQIPJPLRNSP RQIPL RSNRQ RPJQLSNSPRQPRRSSTSVQXPZ RRSSV RPXQ[ RSTRVPXPZQ[S\\"; 
			case '}': return " 38MXQHSITJTLRNQP RSITL RQNRQ RTJSLQNQPRQTRRSQTQVSXTZ RRSQV RTXS[ R QTRVTXTZS[Q\\";
			default: return null;
			}
		};
	
	
Hershey.prototype.charToPathOp=function(letter)
		{
		var i;
		if(letter==' ') return [];
		var s=this.charToHersheyString(letter);
		
		if(s==null) return [];
		
		var num_vertices=0;
		for( i=0;i< 3;++i)
			{
			var c=s.charAt(i);
			if(c==' ' || c=='\n') continue;
			num_vertices = num_vertices*10+(c-'0');
			}
		num_vertices--;
		i+=2;
		var nop=0;
		array=[];
		while(nop<num_vertices)
			{
			var pathOp={operator:null,x:0,y:0};
			pathOp.operator=(array.length==0?Hershey.MOVETO:Hershey.LINETO);
			var c=s.charAt(i++);
			if(c==' ')
				{
				c=s.charAt(i++);
				if(c!='R') throw "BOUM "+s;
				nop++;
				pathOp.operator=Hershey.MOVETO;
				c=s.charAt(i++);
				}
			pathOp.x=c.charCodeAt(0)-'R'.charCodeAt(0);
			c=s.charAt(i++);
			pathOp.y=c.charCodeAt(0)-'R'.charCodeAt(0);
			nop++;
			array.push(pathOp);
			}
		return array;
		};


Hershey.prototype.paint=function(ctx,s,x,y,width,height)
	{
	if(s==null || s.length==0 || width==0 || height==0) return "";
	
	var i,dx=width/s.length;
	for(i=0;i < s.length;++i)
		{
		var n,array=this.charToPathOp(s.charAt(i));
		for(n=0;n< array.length;++n)
			{
			var p2= array[n];
			var x2= x+ (p2.x/this.scalex)*dx + dx*i +dx/2.0;
			var y2= y+ (p2.y/this.scaley)*height +height/2.0 ;
			
			
			
			if(p2.operator == Hershey.LINETO)
				{
				ctx.lineTo(x2,y2);
				}
			else
				{
				ctx.moveTo(x2,y2);
				}
			}
		
		}
	
	};

