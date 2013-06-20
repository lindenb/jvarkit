/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008
Whitehead Institute for Biomedical Research, Steve Rozen
(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky
All rights reserved.

This file was copied from primer3: http://sourceforge.net/p/primer3/code/HEAD/tree/primer3/trunk/src/dpal.c

And transformed to java (although, I didn't test all the dpal )

Used here with the permissions of Dr Rozen and Dr Skaletsky.


    This file is part of primer3 software suite.

    This software suite is is free software;
    you can redistribute it and/or modify it under the terms
    of the GNU General Public License as published by the Free
    Software Foundation; either version 2 of the License, or (at
    your option) any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this software (file gpl-2.0.txt in the source
    distribution); if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
package com.github.lindenb.jvarkit.util.align;

import java.util.Arrays;
import java.util.logging.Logger;


public class Dpal
	{
	int DPAL_MAX_ALIGN=Integer.MAX_VALUE;
	boolean DPAL_FORGET_PATH=false;
	boolean DPAL_PRINT_COVERAGE=false;
	public interface Sequence 
		{
		public int size();
		public char get(int index);
		}
	
	public static class DefaultSequence implements Sequence
		{
		private String s;
		DefaultSequence(String s) { this.s=s;}
		@Override
		public char get(int index) {
			return s.charAt(index);
		}
		@Override
		public int size() {
			return s.length();
		}
		}		
		
	private static final Logger LOG=Logger.getLogger("dpal");
	private static class ExtensibleList
		{
		private int _size=0;
		private int L[]=new int[502];
		public ExtensibleList()
			{
			
			}
		
		public int get(int i)
			{
			return L[i];
			}
		public void set(int index,int val)
			{
			if(index>=L.length)
				{
				LOG.info("extend to "+index);
				int L2[]=new int[index+index/2+1];
				System.arraycopy(L,0,L2, 0,L.length);
				L=L2;
				}
			if(index>=this._size) this._size=index+1;
			L[index]=val;
			}
		}
	
	private static class ExtensibleString
		{
		private int _size=0;
		private char L[]=new char[0];
		public ExtensibleString()
			{
			
			}
		
		public char get(int i)
			{
			return i>=L.length?' ':L[i];
			}
		public void set(int index,char val)
			{
			if(index>=L.length)
				{
				LOG.info("extend to "+index);
				char L2[]=new char[index+index/2+1];
				Arrays.fill(L2, ' ');
				System.arraycopy(L,0,L2, 0,L.length);
				L=L2;
				}
			if(index>=this._size) this._size=index+1;
			L[index]=val;
			}
		}
	
   private static final int DPAL_ERROR_SCORE=Integer.MIN_VALUE;
	  /* 0 means do not exit on error. */

/* 
				       * The maximum size of a string that can be
				       * aligned with with generic dpal and for which
				       * we can return a "path".  Several arrays of
				       * size DPAL_MAX_ALIGN X DPAL_MAX_ALIGN are
				       * statically allocated in dpal.o */

	private static final int DPAL_LOCAL=0;  /* Return a local alignment. */
	private static final int DPAL_GLOBAL_END   =1;  /* 
				      * Return a global alignment _anchored at the end
				      * of the first sequence_.
				      */
	private static final int  DPAL_GLOBAL       =2;  /* 
				      * Return an arbitrary global alignment, that is
				      * one anchored at the end of either the first or
				      * the second sequence.
				      */
	private static final int DPAL_LOCAL_END    =3;   /* 
	                               * Return a local alignment that includes the
				       * end (but not necessarily the beginning) of
				       * the first sequence.
				       */

	/*
	 * It is not possible to specify end-gap penalties for the DPAL_GLOBAL_END
	 * and DPLAL_GLOBAL flags.
	 */

	/* 
	 * The data structure that stores the "scoring system matrix". (The socring
	 * system matrix data structure is of size UCHAR_MAX + 1 by UCHAR_MAX + 1.)
	 */
	private static class dpal_ssm
		{
		public int f(char i,char j)
			{
			if (('A' == i || 'C' == i || 'G' == i || 'T' == i || 'N' == i)
			          && ('A' == j || 'C' == j || 'G' == j || 'T' == j 
			              || 'N' == j)) {
			        if (i == 'N' || j == 'N') 
			          return -25;
			        else if (i == j)
			          return  100;
			        else 
			         return -100;
			      } else
			        return Integer.MIN_VALUE;
			
			}
		}

	/* Structure for passing in arguments to the main function, dpal. */
	static class dpal_args {
	    int check_chars=0;        /* 
			             * If non-0, check for and raise an error on an
			             * illegal character in the input strings.
				     */
	    boolean debug=true;              /* 
				     * If non-0, print debugging information to
				     * stderr.
				     */
	    int fail_stop=1;           /* Exit with -1 on error. */
	    int flag=DPAL_LOCAL;                /* 
				      * One of DPAL_GLOBAL, DPAL_LOCAL,
				      * DPAL_GLOBAL_END, DPAL_LOCAL_END
				      */
	    int force_generic=0;      /* Force the use of the generic function. */
	    int force_long_generic=0; /* 
				     * Force the use of the long generic no-path
				     * function.
				     */
	    int force_long_maxgap1=0; /* Force the use of the long maxgap 1 functions. */
	    int gap=-200;                 /* The "gap opening" penalty. */
	    int gapl=-200;                /* The "gap extension" penalty. */
	    int max_gap=1;             /* 
			              * The maximum allowable size for a gap. -1
			              * indicates that the gap can be of any size.
				      */
	    int score_max;           /* If greater than 0 stop search as soon as
				      * score > score_max.
				      */
	    int score_only=0;          /* 
				      * If non-0, only print the score on
				      * stdout. (Incompatible with debug.)
				      */
	    dpal_ssm ssm=new dpal_ssm();            /* The scoring system matrix. */
	} ;

	/*
	public static class Alignment
		{
		private Sequence seq1;
		private Sequence seq2;
		private int size=0;
		private int array[];
		
		public Sequence getSequence(int y)
			{
			switch(y)
				{
				case 0: return seq1;
				case 1: return seq2;
				default:throw new ArrayIndexOutOfBoundsException("no sequence idx:"+y);
				}
			}		
		
		public void set(int y,int index,int v)
			{
			if(index*2 >= array.length)
				{
				int L2[]=new int[index+index/2+1];
				System.arraycopy(array,0,L2, 0,array.length);
				array=L2;	
				}
			if(index>=this.size) this.size=index+1;
			array[index*2+y]=v;
			}
		
		public int get(int y, int pos)
			{
			return array[pos*2+y];
			}
		
		public int size()
			{
			return this.size;
			}
		}*/
	
	/* Structure for receiving results from the main function, dpal. */
	static class dpal_results {
	    String msg;
	    Matrix2     path=new Matrix2();
	    int     path_length;
	    int     align_end_1; /* Last alignment position in the 1st sequence. */
	    int     align_end_2; /* Last alignment position in the 2nd sequence. */
	    double  score;
	} 


	/* 
	   The next two macros require that the output argument is always
	   called 'out'.
	 */

	static void CHECK_ERROR(boolean COND,String MSG)
		{
		if (COND) { throw new IllegalArgumentException(MSG); }
		}

	/* If fail_stop is set, the code at FAIL: will use fprintf, otherwise
	   will return with msg set. */
	static void DPAL_OOM_ERROR()
		{
		throw new IllegalArgumentException("Out of memory");
		}

	static void fail_action( dpal_args in, dpal_results out) {
	  if (in.fail_stop==1) {
	    System.err.printf( "\n%s\n", out.msg);
	    System.exit(-1);
	  } 
	}

	private  static class Matrix2
		{
		private int array[]=new int[0];
		int num_x=0;
		int num_y=0;
		
		
		private int index(int x,int y)
			{	
			return num_x*y+x;
			}
		int get(int x,int y)
			{
			return array[index(x,y)];
			}
		void set(int x,int y,int v)
			{
			if(x>=num_x || y>=num_y)
				{
				
				int new_x=Math.max(num_x, x+1);
				int new_y=Math.max(num_y, y+1);
				//LOG.info("resize to "+new_x+"/"+new_y);
				int a2[]=new int[new_x*new_y];
				for(int i=0;i< num_x;++i)
					{
					for(int j=0;j< num_y;++j)
						{
						a2[new_x*j+i]=this.get(i,j);
						}
					}
				this.array=a2;
				this.num_x=new_x;
				this.num_y=new_y;
				}
			array[index(x,y)]=v;
			}
		void shift1(int mgy)
			{
			for(int y=0;y< this.num_y;++y) 
				{
				int save=this.get(0,y);
				for(int x=0;x< mgy;++x)
					{
					this.set(x, y,this.get(x+1,y));
					}
				this.set(0,mgy, save);
				}
			}
		}

private  static class Matrix3
	{
	private int array[]=new int[0];
	private int num_x=0;
	private int num_y=0;
	private int num_z=0;

	private int index(int x,int y,int z)
		{	
		int n= num_x*y+x +z*(num_x*num_y);
		if(n<0) throw new ArrayIndexOutOfBoundsException(""+x+"/"+y+" : "+num_x+"/"+num_y);
		return n;
		}

	
	public int get(int x,int y,int z)
			{
			return array[index(x,y,z)];
			}

	public	void set(int x,int y,int z,int v)
			{
			if(x>=num_x || y>=num_y || z>=num_z)
				{
				
				int new_x=Math.max(num_x, x+1);
				int new_y=Math.max(num_y, y+1);
				int new_z=Math.max(num_z, z+1);
				//LOG.info("resize to "+new_x+"/"+new_y+"/"+new_z);
				int a2[]=new int[new_x*new_y*new_z];
				for(int i=0;i< num_x;++i)
					{
					for(int j=0;j< num_y;++j)
						{
						for(int k=0;k< num_z;++k)
							{
							a2[new_x*j+i +k*(new_x*new_y)]=this.get(i,j,k);
							}
						}
					}
				this.array=a2;
				this.num_x=new_x;
				this.num_y=new_y;
				this.num_z=new_z;
				}
			array[index(x,y,z)]=v;

			}
	

	}


void
dpal( /* unsigned */ Sequence X,
      /* unsigned */ Sequence Y,
      dpal_args in,
     dpal_results out)
{
  int xlen, ylen;

  out.score = DPAL_ERROR_SCORE;
  out.path_length = 0;
  out.msg = null;


  xlen =X.size();
  ylen = Y.size();

  out.align_end_1 = -1;
  out.align_end_2 = -1;

  if (xlen==0) {
    out.msg = "Empty first sequence";
    out.score = 0;
    return;
  }
  if (ylen==0) {
    out.msg = "Empty second sequence";
    out.score = 0;
    return;
  }
  if (1 == in.force_generic || in.debug  || 0 == in.score_only) {
    /* 
     * A true value of in->debug really means "print alignment on stderr"
     * and implies 0 == a.score_only.
     */
	  LOG.info("iCI");
    _dpal_generic(X, Y, in, out);
  } else if (1 == in.force_long_generic) {
    _dpal_long_nopath_generic(X, Y, in, out);
  } else if (1 == in.max_gap ) {
    if (DPAL_LOCAL == in.flag)
      _dpal_long_nopath_maxgap1_local(X, Y, in, out);
    else if (DPAL_GLOBAL_END == in.flag)
      _dpal_long_nopath_maxgap1_global_end(X, Y,  in, out);
    else if (DPAL_LOCAL_END == in.flag)
      _dpal_long_nopath_maxgap1_local_end(X, Y, in, out);
    else if (xlen <= DPAL_MAX_ALIGN && ylen <= DPAL_MAX_ALIGN)
      _dpal_generic(X, Y,  in, out);
    else _dpal_long_nopath_generic(X, Y, in, out);
  }
  else if (xlen < DPAL_MAX_ALIGN && ylen < DPAL_MAX_ALIGN)
    _dpal_generic(X, Y,  in, out);
  else
    _dpal_long_nopath_generic(X, Y, in, out);
}

 void
_dpal_generic( /* unsigned */ Sequence X,
               /* unsigned */ Sequence Y,
               
               dpal_args in,
              dpal_results out)
{
	int xlen=X.size(); 
    int ylen=Y.size();
    Matrix2 S=new Matrix2();
   Matrix3 P=null;
if(!DPAL_FORGET_PATH){
  P=new Matrix3();
}

     int i, j, k=0, mg, c;
     int gap = in.gap, gapl = in.gapl, max_gap = in.max_gap;


    int i0 = -99, j0 = -99;
    int saved_k=0;


    int I = -99, J = -99; /* Coordinates of the maximum score. */
    int smax;             /* The optimum score. */
    int score = -99;      /* Current score. */

    int a,b,max;

LOG.info("_dpal_generic called\n");



    /* Initialize the 0th column of the score matrix. */
    smax = Integer.MIN_VALUE;
    for(i=0; i < xlen; i++) {
        score = in.ssm.f(X.get(i),Y.get(0)); 
        if (DPAL_LOCAL == in.flag) {
            if (score < 0) score = 0;
            if(score > smax) {
                smax = score;
                I=i; J=0;
            }
        }
        else if (DPAL_LOCAL_END == in.flag) {if (score < 0) score = 0;}
        S.set(i,0,score);
    }   
    /* Move code for find global-alignment and end-anchored
       alignment below? */
    if (DPAL_LOCAL != in.flag) {
        /* 
         * For a non-local alignment we restrict our search for the maximum
         * score to the last row.
         */
        smax = S.get(xlen-1,0); I=xlen-1; J=0;
    }
           
    /* Initialize the 0th row of the score matrix. */
    for(j=0; j<ylen; j++) { 
        score = in.ssm.f(X.get(0),Y.get(j)); 

        if(DPAL_LOCAL == in.flag){
            if (score < 0) score = 0;
            if(score > smax){
                smax = score;
                I=0; J=j;
            }
        }
        else if (DPAL_LOCAL_END == in.flag) {if (score < 0) score = 0;}
        S.set(0,j,score);
    }   
    if(DPAL_GLOBAL == in.flag&&S.get(0,ylen-1)>smax){
                smax = S.get(0,ylen-1);
                I=0; J=ylen-1;
    }

    /* Further is the solution for dynamic programming problem. */
    for(i=1; i<xlen; i++) {
        for(j=1; j<ylen; j++) {

            a=S.get(i-1,j-1);

            b = c = Integer.MIN_VALUE;
            if (1 == max_gap) {
                if (i > 1) {
                    b = S.get(i-2,j-1) + gap;
if(!DPAL_FORGET_PATH){
                    i0 = i - 2;
}
                }
                if (j > 1) {
                    c = S.get(i-1,j-2) + gap;
if(!DPAL_FORGET_PATH){
                    j0 = j - 2;
}
                }
            } else if (max_gap > 1) {
                max=Integer.MIN_VALUE;
                mg=(max_gap+1>i||max_gap<0)?i:max_gap+1;
                for(k=2; k<=mg; k++) {
                    c = S.get(i-k,j-1) + gap + gapl*(k-2);
                    if(c>max){
                        max=c;
if(!DPAL_FORGET_PATH){
                        i0 = i-k;
}
                    }
                }
                b=max;

                max=Integer.MIN_VALUE;
                mg=(max_gap+1>j||max_gap<0)?j:max_gap+1;
                for(k=2;k<=mg;k++) {
                    c = S.get(i-1,j-k) + gap + gapl*(k-2);
                    if(c>max){
                        max=c;
if(!DPAL_FORGET_PATH){
                        j0 = j-k;
}
                    }
                }
                c=max;
            }

            if(a>=b && a>=c) {
                score = a + in.ssm.f(X.get(i),Y.get(j));
if(!DPAL_FORGET_PATH){
                P.set(i,j,1,i-1);
                P.set(i,j,2,j-1);
}
            } else if (b > a && b >= c) {
                score = b + in.ssm.f(X.get(i),Y.get(j));
if(!DPAL_FORGET_PATH){
                P.set(i,j,1,i0);
                P.set(i,j,2,j-1);
}
            } else if (c > a && c > b) {
                score = c + in.ssm.f(X.get(i),Y.get(j));
if(!DPAL_FORGET_PATH){
                P.set(i,j,1,i-1);
                P.set(i,j,2,j0);
}
            }

            if (score >= smax)
                /* 
                 * Because of comparison '>=' immediately above, dpal reports
                 * ungapped (i.e. diagonal) alignments if there is a choice
                 * of more than one optimum alignment.
                 */
                /* Move code to get 'g' and 'e' maxima to a separate loop ? */
                if (DPAL_LOCAL == in.flag 
                    || (DPAL_GLOBAL_END == in.flag && i == xlen-1)
                    || (DPAL_LOCAL_END   == in.flag && i == xlen-1)
                    || (DPAL_GLOBAL == in.flag&& (i==xlen-1||j==ylen-1))) {
                    /*  
                     * If in->flag is DPAL_LOCAL, then a cell anywhere within
                     * S may be the endpoint of the alignment.  If in->flag is
                     * DPAL_GLOBAL_END, then only cells in the last row may be
                     * the endpoint.  If in->flag is DPAL_GLOBAL cells in the
                     * last row _or_ the last column may be the endpoint.
                     */
                    smax = score;
                    I = i;
                    J = j;
            } /*  put else here ? */
            if (score < 0 && (DPAL_LOCAL == in.flag 
                              || DPAL_LOCAL_END == in.flag))
                /* 
                 * For a local alignment, 0 is the lowest score that we record
                 * in S.
                 */
                score = 0;

            S.set(i,j,score);
        }
    }
    /* I and J now specify the last pair of an optimum alignment. */

if(!DPAL_FORGET_PATH){    
    k = (I > J) ? I+1 : J+1;
    saved_k=k;

    out.path.set(k,0,I); out.path.set(k,1,J);
    while(out.path.get(k,0)!=0&&out.path.get(k,1)!=0) {
        if ((in.flag== DPAL_LOCAL || in.flag == DPAL_LOCAL_END)
                 &&S.get(out.path.get(k,0),out.path.get(k,1))==0) {
          k++; break;
        }
        out.path.set(k-1,0,P.get(out.path.get(k,0),out.path.get(k,1),1));
        out.path.set(k-1,1,P.get(out.path.get(k,0),out.path.get(k,1),2));
        k--;
    }
    if (k>0) {
        for (i=0;i<=saved_k-k;i++) {
            out.path.set(i,0,out.path.get(i+k,0));
            out.path.set(i,1,out.path.get(i+k,1));
        }
    }
}

    if ((DPAL_LOCAL == in.flag 
         || DPAL_LOCAL_END == in.flag)&& S.get(I,J) <= 0) {
        /* There is no alignment at all. */
        out.score = 0;
        out.path_length = 0;
    } else {
        out.score = smax;
        out.align_end_1 = I;
        out.align_end_2 = J;
if(!DPAL_FORGET_PATH){        
        out.path_length = saved_k - k + 1;
}
else
{
        out.path_length = 0;
}
    }
    LOG.info("Donexx");
if(!DPAL_FORGET_PATH){LOG.info("Done");
    if (in.debug) print_align(X,Y,P,I,J, in);
}
    return;

} /* _dpal_generic */

/* Linear space, no path, for any value of maxgap and for any alignment. */
 void
_dpal_long_nopath_generic( /* unsigned */ Sequence X,
                           /* unsigned */ Sequence Y,
                          
                           dpal_args in,
                          dpal_results out)
{
	
	 int xlen=X.size();
     int ylen=Y.size();

    /* The "score matrix" (matrix of best scores). */
    Matrix2 S=new Matrix2();
    //ExtensibleList SI=new ExtensibleList();
    //Matrix2 P=new Matrix2();    

     int i, j, k, mg, mgy, c;
     int gap = in.gap, gapl = in.gapl, max_gap = in.max_gap;


    int I = -99, J = -99; /* Coordinates of the maximum score. */
    int smax;             /* The optimum score. */
    int score;            /* Current score. */


    out.score = DPAL_ERROR_SCORE;
    out.path_length = 0;
    out.msg = null;


    /* Initialize the 0th column of the score matrix. */
    smax = Integer.MIN_VALUE;
    for(i=0; i < xlen; i++) {
        score = in.ssm.f(X.get(i),Y.get(0)); 
        if (DPAL_LOCAL == in.flag) {
            if (score < 0) score = 0;
            if(score > smax) {
                smax = score;
                I=i; J=0;
            }
        }
        else if (DPAL_LOCAL_END == in.flag) {if (score < 0) score = 0;}
        S.set(0,i,score);
    }   
    /* Move code for find global-alignment and end-anchored
       alignment below? */
    if (DPAL_LOCAL != in.flag) {
        /* 
         * For a non-local alignment we restrict our search for the maximum
         * score to the last row.
         */
        smax = S.get(0,xlen-1); I=xlen-1; J=0;
    }
           
    /* Initialize the 0th row of the score matrix. 
    for(j=0; j<ylen; j++) { 
        score = in->ssm[X[0]][Y[j]]; 

        if(DPAL_LOCAL == in->flag){
            if (score < 0) score = 0;
            if(score > smax){
                smax = score;
                I=0; J=j;
            }
        }
        S[0][j] = score;
    }   
    
    if(DPAL_GLOBAL == in->flag&&S[0][ylen-1]>smax){
                smax = S[0][ylen-1];
                I=0; J=ylen-1;
    }
    */

    /* Further is the solution for dynamic programming problem. */
    for(j=1; j<ylen; j++) {
        mgy=(max_gap+1>j||max_gap<0)?j:max_gap+1;
        score = in.ssm.f(X.get(0),Y.get(j));
         if (DPAL_LOCAL == in.flag) {
             if (score < 0) score = 0;
             if(score > smax) smax = score;
         }    
         else if (DPAL_LOCAL_END == in.flag) { if (score < 0) score = 0;}
         else if (DPAL_GLOBAL == in.flag && j == ylen-1 && score > smax)
                        smax = score;
        S.set(mgy,0,score);
        for(i=1; i<xlen; i++) {

            score=S.get(mgy-1,i-1);

                mg=(max_gap+1>i||max_gap<0)?i:max_gap+1;
                for(k=2; k<=mg; k++) 
                    if((c = S.get(mgy-1,i-k) + gap + gapl*(k-2)) > score)score = c;

                for(k=2;k<=mgy;k++) 
                    if((c = S.get(mgy-k,i-1) + gap + gapl*(k-2)) > score)score=c;

                score += in.ssm.f(X.get(i),Y.get(j));

            if (score >= smax)
                /* 
                 * Because of comparison '>=' immediately above, dpal reports
                 * ungapped (i.e. diagonal) alignments if there is a choice
                 * of more than one optimum alignment.
                 */
                /* Move code to get 'g' and 'e' maxima to a separate loop ? */
                if (DPAL_LOCAL == in.flag 
                    || ((DPAL_GLOBAL_END == in.flag
                         || DPAL_LOCAL_END == in.flag) 
                        && i == xlen-1)
                    || (DPAL_GLOBAL == in.flag&& (i==xlen-1||j==ylen-1))) {
                    /*  
                     * If in->flag is DPAL_LOCAL, then a cell anywhere within
                     * S may be the endpoint of the alignment.  If in->flag is
                     * DPAL_GLOBAL_END, then only cells in the last row may be
                     * the endpoint.  If in->flag is DPAL_GLOBAL cells in the
                     * last row _or_ the last column may be the endpoint.
                     */
                    smax = score;
                    I = i;
                    J = j;
            } /*  put else here ? */
            if (score < 0 && (DPAL_LOCAL == in.flag
                              || DPAL_LOCAL_END == in.flag))
                /* 
                 * For a local alignment, 0 is the lowest score that we record
                 * in S.
                 */
                score = 0;

            S.set(mgy,i,score);
        }
        if(mgy == max_gap + 1){
        	S.shift1(mgy);
        }
    }
    /* I and J now specify the last pair of an optimum alignment. */

    if (DPAL_LOCAL == in.flag && smax <= 0) {
        /* There is no alignment at all. */
        out.score = 0;
        out.path_length = 0;
    } else {
        out.score = smax;
        out.align_end_1 = I;
        out.align_end_2 = J;
    }
} /* _dpal_long_nopath_generic */

 void
_dpal_long_nopath_maxgap1_local( /* unsigned */ Sequence X,
                                 /* unsigned */ Sequence Y,
                                
                                 dpal_args in,
                                dpal_results out)
{
	 int xlen=X.size(); 
     int ylen=Y.size();
     ExtensibleList S0, S1, S2; 
    ExtensibleList P0=new ExtensibleList();
    ExtensibleList P1=new ExtensibleList();
    ExtensibleList P2=new ExtensibleList();
    
    
    ExtensibleList S=new ExtensibleList();

     int i, j;
     int gap = in.gap;
     int smax;           /* The optimum score. */
     int score;          /* Current score. */
     int a;




    S0 = P0; S1 = P1; S2 = P2;

    smax = 0; /* For local alignment score can never be less than 0. */

    /* Initialize the 0th row of the score matrix. */
    for(j=0; j < ylen; j++) { 
        score = in.ssm.f(X.get(0),Y.get(j)); 
        if (score < 0) score = 0;
        else if (score > smax) smax = score;
        /*S[0][j] = score;*/
        S0.set(j,score);
    }   

    /* Set the 1st row of the score matrix. */
    score = in.ssm.f(X.get(1),Y.get(0));
    if(score < 0) score = 0;
    else if (score > smax) smax = score;
    S1.set(0,score);
    for(j=1; j < ylen; j++) {
        score = S0.get(j-1);
        if(j>1 && (a=S0.get(j-2) + gap) > score)score = a;
        score += in.ssm.f(X.get(1),Y.get(j));
        if (score < 0) score = 0;
        else if(score > smax) smax = score;
        S1.set(j,score);
    }

    for(i=2; i < xlen; i++) {
        score = in.ssm.f(X.get(i),Y.get(0));
        if (score < 0) score = 0;
        else if (score > smax) smax = score;
        S2.set(0,score);
        score = S1.get(0);
        if((a=S0.get(0) + gap) > score) score = a;
        score += in.ssm.f(X.get(i),Y.get(1));
        if(score < 0) score = 0;
        else if (score > smax) smax = score;
        S2.set(1,score);
        for(j=2; j < ylen; j++) {
            score = S0.get(j-1);
            if((a=S1.get(j-2))>score) score = a;
            score +=gap;
            if((a=S1.get(j-1)) >score) score = a;

            score += in.ssm.f(X.get(i),Y.get(j));       
            if (score < 0 ) score = 0;
            else if (score > smax) smax = score;
            S2.set(j,score);
        }
        S = S0; S0 = S1; S1 = S2; S2 = S;
    }
    out.score = smax;
    out.path_length=0;
    P0=null; P1=null; P2=null;
    return;
} /* _dpal_long_nopath_maxgap1_local */

void
_dpal_long_nopath_maxgap1_global_end( /* unsigned */ Sequence X,
        /* unsigned */ Sequence Y,
                                    
                                      dpal_args in,
                                     dpal_results out)
{
	  int xlen=X.size(); 
      int ylen=Y.size();
  /* The "score matrix" (matrix of best scores). */
      ExtensibleList P0=new ExtensibleList();
      ExtensibleList P1=new ExtensibleList();
      ExtensibleList P2=new ExtensibleList();
      ExtensibleList S,S0,S1,S2;

   int i, j, k;
   int gap = in.gap;
   int smax;           /* The optimum score. */
   int score;          /* Current score. */
   int a, t;


  S0 = P0; S1 = P1; S2 = P2;

  smax = in.ssm.f(X.get(xlen-1),Y.get(0));
           
  /* Set the 0th row of the score matrix. */
  for(j=0; j<xlen; j++) S0.set(j,in.ssm.f(X.get(j),Y.get(0)));

  /* Set the 1st row of the score matrix. */
  S1.set(0,in.ssm.f(X.get(0),Y.get(1)));
  for(j=1; j < xlen; j++){
    score = S0.get(j-1);
    if(j>1 && (a=S0.get(j-2) + gap)> score)score = a;
    score += in.ssm.f(X.get(j),Y.get(1));
    if(score > smax && j == xlen-1) smax = score;
    S1.set(j,score);
  }

  k = ylen - (int)(xlen / 2) + 1;
  if (k<1) k = 1;

  /* Set the rectangular part of almost the remainder of the matrix. */
  for(j=2; j<k+1; j++) {
    S2.set(0,in.ssm.f(X.get(0),Y.get(j)));
    score = S1.get(0);
    if((a=S0.get(0)+gap) > score) score = a;
    score += in.ssm.f(X.get(1),Y.get(j));
    S2.set(1,score);
    for(i=2; i<xlen-1; i++) {
      score = S1.get(i-2);
      if((a=S0.get(i-1)) > score)score = a;
      score += gap;
      if((a=S1.get(i-1)) > score)score = a;
      score += in.ssm.f(X.get(i),Y.get(j));
      S2.set(i,score);
    }
    score = S1.get(xlen-3);
    if((a=S0.get(xlen-2)) > score)score = a;
    score += gap;
    if((a=S1.get(xlen-2)) > score)score = a;
    score += in.ssm.f(X.get(xlen-1),Y.get(j));
    S2.set(xlen-1,score);
    if(score > smax) smax = score;
    S = S0; S0 = S1; S1 = S2; S2 = S;
  }

  /* Set the triangular part of almost the remainder of the matrix. */
  t = 2;
  for(j=k+1; j<ylen; j++) {
    for(i=t; i<xlen-1; i++) {
      score = S1.get(i-2);
      if((a=S0.get(i-1)) > score) score = a;
      score += gap;
      if((a=S1.get(i-1)) > score) score = a;
      score += in.ssm.f(X.get(i),Y.get(j));
      S2.set(i,score);
    }
    t += 2;
    score = S1.get(xlen-3);
    if((a=S0.get(xlen-2)) > score)score = a;
    score += gap;
    if((a=S1.get(xlen-2)) > score)score = a;
    score += in.ssm.f(X.get(xlen-1),Y.get(j));
    S2.set(xlen-1,score);
    if(score > smax) smax = score;
    S = S0; S0 = S1; S1 = S2; S2 = S;
  }

  P0=null; P1=null; P2=null;
  out.score = smax;
  out.path_length=0;
  return;

} /* _dpal_long_nopath_maxgap_global_end */





void print_align( Sequence X,
		Sequence Y,
            Matrix3 P,
            int I,
            int J,
             dpal_args dargs)
{
	ExtensibleList JX=new ExtensibleList();
	ExtensibleList JY=new ExtensibleList();
	ExtensibleString sx=new ExtensibleString();
	ExtensibleString sy=new ExtensibleString();
	ExtensibleString sxy=new ExtensibleString();
	
  int k,i,j,n,m;

  if(I>J)k=I+1;
  else k=J+1;

  n=k;
  JX.set(k,I);
  JY.set(k,J);
  while(JX.get(k)!=0&&JY.get(k)!=0){
    JX.set(k-1,P.get(JX.get(k),JY.get(k),1));
    JY.set(k-1,P.get(JX.get(k),JY.get(k),2));
    k--;
  }
  if(JX.get(k)>JY.get(k)){
    for(i=0;i<JX.get(k);i++)sx.set(i,X.get(i));
    for(i=0;i<JX.get(k)-JY.get(k);i++)sy.set(i,' ');
    j = JX.get(k)-JY.get(k);
    for(i=JX.get(k)-JY.get(k);i<JX.get(k);i++)sy.set(i,Y.get(i-j));
    m = JX.get(k);
  }
  else{
    for(i=0;i<JY.get(k);i++)sy.set(i,Y.get(i));
    for(i=0;i<JY.get(k)-JX.get(k);i++)sx.set(i,' ');
    j= JY.get(k)-JX.get(k);
    for(i=j;i<JY.get(k);i++)sx.set(i,X.get(i-j));
    m = JY.get(k);
  }
  for(i=0;i<m;i++)sxy.set(i,' ');
  for(i=k;i<n;i++){
    sx.set(m,X.get(JX.get(i)));
    sy.set(m,Y.get(JY.get(i)));
    /* if(sx[m]==sy[m]&&sx[m]!='N') sxy[m] = '|'; */
    if (dargs.ssm.f((/* unsigned */ char)sx.get(m),(/* unsigned */ char)sy.get(m)) > 0)
      sxy.set(m,'|');
    else sxy.set(m,' ');
    if(JX.get(i+1)-JX.get(i)>JY.get(i+1)-JY.get(i)){
      for(j=1;j<JX.get(i+1)-JX.get(i);j++){
        sy.set(m+j,'-');
        sx.set(m+j,X.get(JX.get(i)+j));
        sxy.set(m+j,' ');
      }
      m += JX.get(i+1)-JX.get(i)-1;
    }
    if(JY.get(i+1)-JY.get(i)>JX.get(i+1)-JX.get(i)){
      for(j=1;j<JY.get(i+1)-JY.get(i);j++){
        sx.set(m+j,'-');
        sy.set(m+j,Y.get(JY.get(i)+j));
        sxy.set(m+j,' ');
      }
      m += JY.get(i+1)-JY.get(i)-1;
    }
    m++;
  }
  sx.set(m,X.get(I));
  sy.set(m,Y.get(J));
  for (i=m+1; i < (m + X.size() - I); i++) 
    sx.set(i,X.get(i-m+I));
  for (i=m+1; i < (m + Y.size() - J); i++) 
    sy.set(i,Y.get(i-m+J));

  if (dargs.ssm.f((/* unsigned */ char)sx.get(m),(/* unsigned */ char)sy.get(m)) > 0)
    sxy.set(m,'|');
  else sxy.set(m,' ');
  m++;
  if (X.size() - I >Y.size()-J) {
    k = m + X.size() - I;
  } else {
    k = m + Y.size() - J;
  }
        
  j=0;
  while(j<k){
    for(i=j;i<j+70;i++) System.err.printf( "%c",sx.get(i));
    System.err.printf( "\n");
    for(i=j;i<j+70;i++) System.err.printf( "%c",sxy.get(i));
    System.err.printf( "\n");
    for(i=j;i<j+70;i++) System.err.printf( "%c",sy.get(i)); System.err.printf("\n");
    for(i=0;i<70;i++)   System.err.printf( "_");
    System.err.printf( "\n");
    j +=70;
  }
}  /* print_align(X,Y,P,I,J, dargs) */



	
 void
_dpal_long_nopath_maxgap1_local_end( Sequence X,
									Sequence Y,
                                     dpal_args in,
                                    dpal_results out)
{
  /* The "score matrix" (matrix of best scores). */
  ExtensibleList P0=new ExtensibleList(); 
  ExtensibleList P1=new ExtensibleList();
  ExtensibleList P2=new ExtensibleList(); 
  ExtensibleList S0=P0;
  ExtensibleList S1=P1;
  ExtensibleList S2=P2;
  ExtensibleList S;
int ylen=Y.size();
int xlen=X.size();
  
   int i, j;
   int gap = in.gap;
   int smax;           /* The optimum score. */
   int score;          /* Current score. */
   int a;

/*#ifdef DPAL_PRINT_COVERAGE
  System.err.printf( "_dpal_long_nopath_maxgap1_local_end called\n");
#endif*/
  /* Note: S2[0] and S2[1] do not get initialized in this case. */


  smax = 0; /* For local alignment score can never be less than 0. */

  /* Initialize the 0th row of the score matrix. */
  for(j=0; j < ylen; j++) { 
    score = in.ssm.f(X.get(0),Y.get(j)); 
    if (score < 0) score = 0;
    /*S[0][j] = score;*/
    S0.set(j,score);
  }   

  /* Set the 1st row of the score matrix. */
  score = in.ssm.f(X.get(1),Y.get(0));
  if(score < 0) score = 0;
  S1.set(0,score);
  for(j=1; j < ylen; j++) {
    score = S0.get(j-1);
    if(j>1 && (a=S0.get(j-2) + gap) > score)score = a;
    score += in.ssm.f(X.get(1),Y.get(j));
    if (score < 0) score = 0;
    S1.set(j,score);
  }

  for(i=2; i < xlen - 1; i++) {
    score = in.ssm.f(X.get(i),Y.get(0));
    if (score < 0) score = 0;
    S2.set(0,score);
    score = S1.get(0);
    if((a=S0.get(0) + gap) > score) score = a;
    score += in.ssm.f(X.get(i),Y.get(1));
    if(score < 0) score = 0;
    S2.set(1,score);
    for(j=2; j < ylen; j++) {
      score = S0.get(j-1);
      if((a=S1.get(j-2))>score) score = a;
      score +=gap;
      if((a=S1.get(j-1)) >score) score = a;

      score += in.ssm.f(X.get(i),Y.get(j));       
      if (score < 0 ) score = 0;
      S2.set(j,score);
    }
    S = S0; S0 = S1; S1 = S2; S2 = S;
  }
  /* Calculate scores for last row (i = xlen-1) and find smax */
  i = xlen - 1;
  score = in.ssm.f(X.get(i),Y.get(0));
  if (score < 0) score = 0;
  else if (score > smax) smax = score;
  S2.set(0,score);
  score = S1.get(0);
  if((a=S0.get(0) + gap) > score) score = a;
  score += in.ssm.f(X.get(i),Y.get(1));
  if(score < 0) score = 0;
  else if (score > smax) smax = score;
  S2.set(1,score);
  for(j=2; j < ylen; j++) {
    score = S0.get(j-1);
    if((a=S1.get(j-2))>score) score = a;
    score +=gap;
    if((a=S1.get(j-1)) >score) score = a;
    score += in.ssm.f(X.get(i),Y.get(j));
    if (score < 0 ) score = 0;
    else if (score > smax) smax = score;
    S2.set(j,score);
  }
  out.score = smax;
  out.path_length=0;
  P0=null; P1=null; P2=null;


} 

	public static void main(String[] args) 
		{
		DefaultSequence X=new DefaultSequence("ACTAGCTGATCGTACGATCGTAACATGCTACGATCGATCGATCA");
		DefaultSequence Y=new DefaultSequence("CTAGCTGATCGTACGATGCTACGATCGATCATCAATGGGCTAGCTGATCGTACG");
		dpal_results r=new dpal_results();
		Dpal app=new Dpal();
		dpal_args a=new dpal_args();
	
		app.dpal(X, Y, a,r);
		System.out.println("Done");
		}


}
