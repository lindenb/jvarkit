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
package com.github.lindenb.jvarkit.util.swing;

import java.awt.Color;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.ParameterException;
import com.github.lindenb.jvarkit.lang.CharSplitter;

import htsjdk.samtools.SAMRecord;


public class ColorUtils
	{

	 /**
	  * YC tag used by the UCSC/IGV
	  * 
	  * See http://software.broadinstitute.org/software/igv/book/export/html/6
	  * http://genome.ucsc.edu/goldenPath/help/hgBamTrackHelp.html
	  */
	public static final String YC_TAG="YC";

	
	public static class Converter implements IStringConverter<Color> {
		public static final String OPT_DESC=" A named color ('red', 'blue'...) use the syntax 'rgb(int,int,int)'.";

		private final ColorUtils cols= new ColorUtils();
		@Override
		public Color convert(final String s) {
			try {
				return this.cols.parse(s);
				}
			catch(Exception err)
				{
				throw new ParameterException("Cannot convert "+s+" to color", err);
				}
			}
		}
	

	public static final Color aliceblue = new Color(240, 248, 255);
	public static final Color antiquewhite = new Color(250, 235, 215);
	public static final Color aquamarine = new Color(127, 255, 212);
	public static final Color aqua = new Color( 0, 255, 255);
	public static final Color azure = new Color(240, 255, 255);
	public static final Color beige = new Color(245, 245, 220);
	public static final Color bisque = new Color(255, 228, 196);
	public static final Color blanchedalmond = new Color(255, 235, 205);
	public static final Color blueviolet = new Color(138, 43, 226);
	public static final Color brown = new Color(165, 42, 42);
	public static final Color burlywood = new Color(222, 184, 135);
	public static final Color cadetblue = new Color( 95, 158, 160);
	public static final Color chartreuse = new Color(127, 255, 0);
	public static final Color chocolate = new Color(210, 105, 30);
	public static final Color coral = new Color(255, 127, 80);
	public static final Color cornflowerblue = new Color(100, 149, 237);
	public static final Color cornsilk = new Color(255, 248, 220);
	public static final Color crimson = new Color(220, 20, 60);
	public static final Color darkblue = new Color( 0, 0, 139);
	public static final Color darkcyan = new Color( 0, 139, 139);
	public static final Color darkgoldenrod = new Color(184, 134, 11);
	public static final Color darkgray = new Color(169, 169, 169);
	public static final Color darkgreen = new Color( 0, 100, 0);
	public static final Color darkgrey = new Color(169, 169, 169);
	public static final Color darkkhaki = new Color(189, 183, 107);
	public static final Color darkmagenta = new Color(139, 0, 139);
	public static final Color darkolivegreen = new Color( 85, 107, 47);
	public static final Color darkorange = new Color(255, 140, 0);
	public static final Color darkorchid = new Color(153, 50, 204);
	public static final Color darkred = new Color(139, 0, 0);
	public static final Color darksalmon = new Color(233, 150, 122);
	public static final Color darkseagreen = new Color(143, 188, 143);
	public static final Color darkslateblue = new Color( 72, 61, 139);
	public static final Color darkslategray = new Color( 47, 79, 79);
	public static final Color darkslategrey = new Color( 47, 79, 79);
	public static final Color darkturquoise = new Color( 0, 206, 209);
	public static final Color darkviolet = new Color(148, 0, 211);
	public static final Color deeppink = new Color(255, 20, 147);
	public static final Color deepskyblue = new Color( 0, 191, 255);
	public static final Color dimgray = new Color(105, 105, 105);
	public static final Color dimgrey = new Color(105, 105, 105);
	public static final Color dodgerblue = new Color( 30, 144, 255);
	public static final Color firebrick = new Color(178, 34, 34);
	public static final Color floralwhite = new Color(255, 250, 240);
	public static final Color forestgreen = new Color( 34, 139, 34);
	public static final Color fuchsia = new Color(255, 0, 255);
	public static final Color gainsboro = new Color(220, 220, 220);
	public static final Color ghostwhite = new Color(248, 248, 255);
	public static final Color goldenrod = new Color(218, 165, 32);
	public static final Color gold = new Color(255, 215, 0);
	public static final Color gray = new Color(128, 128, 128);
	public static final Color greenyellow = new Color(173, 255, 47);
	public static final Color grey = new Color(128, 128, 128);
	public static final Color honeydew = new Color(240, 255, 240);
	public static final Color hotpink = new Color(255, 105, 180);
	public static final Color indianred = new Color(205, 92, 92);
	public static final Color indigo = new Color( 75, 0, 130);
	public static final Color ivory = new Color(255, 255, 240);
	public static final Color khaki = new Color(240, 230, 140);
	public static final Color lavenderblush = new Color(255, 240, 245);
	public static final Color lavender = new Color(230, 230, 250);
	public static final Color lawngreen = new Color(124, 252, 0);
	public static final Color lemonchiffon = new Color(255, 250, 205);
	public static final Color lightblue = new Color(173, 216, 230);
	public static final Color lightcoral = new Color(240, 128, 128);
	public static final Color lightcyan = new Color(224, 255, 255);
	public static final Color lightgoldenrodyellow = new Color(250, 250, 210);
	public static final Color lightgray = new Color(211, 211, 211);
	public static final Color lightgreen = new Color(144, 238, 144);
	public static final Color lightgrey = new Color(211, 211, 211);
	public static final Color lightpink = new Color(255, 182, 193);
	public static final Color lightsalmon = new Color(255, 160, 122);
	public static final Color lightseagreen = new Color( 32, 178, 170);
	public static final Color lightskyblue = new Color(135, 206, 250);
	public static final Color lightslategray = new Color(119, 136, 153);
	public static final Color lightslategrey = new Color(119, 136, 153);
	public static final Color lightsteelblue = new Color(176, 196, 222);
	public static final Color lightyellow = new Color(255, 255, 224);
	public static final Color limegreen = new Color( 50, 205, 50);
	public static final Color lime = new Color( 0, 255, 0);
	public static final Color linen = new Color(250, 240, 230);
	public static final Color magenta = new Color(255, 0, 255);
	public static final Color maroon = new Color(128, 0, 0);
	public static final Color mediumaquamarine = new Color(102, 205, 170);
	public static final Color mediumblue = new Color( 0, 0, 205);
	public static final Color mediumorchid = new Color(186, 85, 211);
	public static final Color mediumpurple = new Color(147, 112, 219);
	public static final Color mediumseagreen = new Color( 60, 179, 113);
	public static final Color mediumslateblue = new Color(123, 104, 238);
	public static final Color mediumspringgreen = new Color( 0, 250, 154);
	public static final Color mediumturquoise = new Color( 72, 209, 204);
	public static final Color mediumvioletred = new Color(199, 21, 133);
	public static final Color midnightblue = new Color( 25, 25, 112);
	public static final Color mintcream = new Color(245, 255, 250);
	public static final Color mistyrose = new Color(255, 228, 225);
	public static final Color moccasin = new Color(255, 228, 181);
	public static final Color navajowhite = new Color(255, 222, 173);
	public static final Color navy = new Color( 0, 0, 128);
	public static final Color oldlace = new Color(253, 245, 230);
	public static final Color olivedrab = new Color(107, 142, 35);
	public static final Color olive = new Color(128, 128, 0);
	public static final Color orangered = new Color(255, 69, 0);
	public static final Color orchid = new Color(218, 112, 214);
	public static final Color palegoldenrod = new Color(238, 232, 170);
	public static final Color palegreen = new Color(152, 251, 152);
	public static final Color paleturquoise = new Color(175, 238, 238);
	public static final Color palevioletred = new Color(219, 112, 147);
	public static final Color papayawhip = new Color(255, 239, 213);
	public static final Color peachpuff = new Color(255, 218, 185);
	public static final Color peru = new Color(205, 133, 63);
	public static final Color plum = new Color(221, 160, 221);
	public static final Color powderblue = new Color(176, 224, 230);
	public static final Color purple = new Color(128, 0, 128);
	public static final Color rosybrown = new Color(188, 143, 143);
	public static final Color royalblue = new Color( 65, 105, 225);
	public static final Color saddlebrown = new Color(139, 69, 19);
	public static final Color salmon = new Color(250, 128, 114);
	public static final Color sandybrown = new Color(244, 164, 96);
	public static final Color seagreen = new Color( 46, 139, 87);
	public static final Color seashell = new Color(255, 245, 238);
	public static final Color sienna = new Color(160, 82, 45);
	public static final Color silver = new Color(192, 192, 192);
	public static final Color skyblue = new Color(135, 206, 235);
	public static final Color slateblue = new Color(106, 90, 205);
	public static final Color slategray = new Color(112, 128, 144);
	public static final Color slategrey = new Color(112, 128, 144);
	public static final Color snow = new Color(255, 250, 250);
	public static final Color springgreen = new Color( 0, 255, 127);
	public static final Color steelblue = new Color( 70, 130, 180);
	public static final Color tan = new Color(210, 180, 140);
	public static final Color teal = new Color( 0, 128, 128);
	public static final Color thistle = new Color(216, 191, 216);
	public static final Color tomato = new Color(255, 99, 71);
	public static final Color turquoise = new Color( 64, 224, 208);
	public static final Color violet = new Color(238, 130, 238);
	public static final Color wheat = new Color(245, 222, 179);
	public static final Color whitesmoke = new Color(245, 245, 245);
	public static final Color yellowgreen = new Color(154, 205, 50);

	
	
    private Map<String,Color> text2color=null;
    public ColorUtils()
		{
        text2color= new HashMap<String,Color>(200);
        text2color.put("aliceblue",aliceblue);
        text2color.put("antiquewhite",antiquewhite);
        text2color.put("aquamarine",aquamarine);
        text2color.put("aqua",aqua);
        text2color.put("azure",azure);
        text2color.put("beige",beige);
        text2color.put("bisque",bisque);
        text2color.put("blanchedalmond",blanchedalmond);
        text2color.put("blueviolet",blueviolet);
        text2color.put("brown",brown);
        text2color.put("burlywood",burlywood);
        text2color.put("cadetblue",cadetblue);
        text2color.put("chartreuse",chartreuse);
        text2color.put("chocolate",chocolate);
        text2color.put("coral",coral);
        text2color.put("cornflowerblue",cornflowerblue);
        text2color.put("cornsilk",cornsilk);
        text2color.put("crimson",crimson);
        text2color.put("darkblue",darkblue);
        text2color.put("darkcyan",darkcyan);
        text2color.put("darkgoldenrod",darkgoldenrod);
        text2color.put("darkgray",darkgray);
        text2color.put("darkgreen",darkgreen);
        text2color.put("darkgrey",darkgrey);
        text2color.put("darkkhaki",darkkhaki);
        text2color.put("darkmagenta",darkmagenta);
        text2color.put("darkolivegreen",darkolivegreen);
        text2color.put("darkorange",darkorange);
        text2color.put("darkorchid",darkorchid);
        text2color.put("darkred",darkred);
        text2color.put("darksalmon",darksalmon);
        text2color.put("darkseagreen",darkseagreen);
        text2color.put("darkslateblue",darkslateblue);
        text2color.put("darkslategray",darkslategray);
        text2color.put("darkslategrey",darkslategrey);
        text2color.put("darkturquoise",darkturquoise);
        text2color.put("darkviolet",darkviolet);
        text2color.put("deeppink",deeppink);
        text2color.put("deepskyblue",deepskyblue);
        text2color.put("dimgray",dimgray);
        text2color.put("dimgrey",dimgrey);
        text2color.put("dodgerblue",dodgerblue);
        text2color.put("firebrick",firebrick);
        text2color.put("floralwhite",floralwhite);
        text2color.put("forestgreen",forestgreen);
        text2color.put("fuchsia",fuchsia);
        text2color.put("gainsboro",gainsboro);
        text2color.put("ghostwhite",ghostwhite);
        text2color.put("goldenrod",goldenrod);
        text2color.put("gold",gold);
        text2color.put("gray",gray);
        text2color.put("greenyellow",greenyellow);
        text2color.put("grey",grey);
        text2color.put("honeydew",honeydew);
        text2color.put("hotpink",hotpink);
        text2color.put("indianred",indianred);
        text2color.put("indigo",indigo);
        text2color.put("ivory",ivory);
        text2color.put("khaki",khaki);
        text2color.put("lavenderblush",lavenderblush);
        text2color.put("lavender",lavender);
        text2color.put("lawngreen",lawngreen);
        text2color.put("lemonchiffon",lemonchiffon);
        text2color.put("lightblue",lightblue);
        text2color.put("lightcoral",lightcoral);
        text2color.put("lightcyan",lightcyan);
        text2color.put("lightgoldenrodyellow",lightgoldenrodyellow);
        text2color.put("lightgray",lightgray);
        text2color.put("lightgreen",lightgreen);
        text2color.put("lightgrey",lightgrey);
        text2color.put("lightpink",lightpink);
        text2color.put("lightsalmon",lightsalmon);
        text2color.put("lightseagreen",lightseagreen);
        text2color.put("lightskyblue",lightskyblue);
        text2color.put("lightslategray",lightslategray);
        text2color.put("lightslategrey",lightslategrey);
        text2color.put("lightsteelblue",lightsteelblue);
        text2color.put("lightyellow",lightyellow);
        text2color.put("limegreen",limegreen);
        text2color.put("lime",lime);
        text2color.put("linen",linen);
        text2color.put("magenta",magenta);
        text2color.put("maroon",maroon);
        text2color.put("mediumaquamarine",mediumaquamarine);
        text2color.put("mediumblue",mediumblue);
        text2color.put("mediumorchid",mediumorchid);
        text2color.put("mediumpurple",mediumpurple);
        text2color.put("mediumseagreen",mediumseagreen);
        text2color.put("mediumslateblue",mediumslateblue);
        text2color.put("mediumspringgreen",mediumspringgreen);
        text2color.put("mediumturquoise",mediumturquoise);
        text2color.put("mediumvioletred",mediumvioletred);
        text2color.put("midnightblue",midnightblue);
        text2color.put("mintcream",mintcream);
        text2color.put("mistyrose",mistyrose);
        text2color.put("moccasin",moccasin);
        text2color.put("navajowhite",navajowhite);
        text2color.put("navy",navy);
        text2color.put("oldlace",oldlace);
        text2color.put("olivedrab",olivedrab);
        text2color.put("olive",olive);
        text2color.put("orangered",orangered);
        text2color.put("orchid",orchid);
        text2color.put("palegoldenrod",palegoldenrod);
        text2color.put("palegreen",palegreen);
        text2color.put("paleturquoise",paleturquoise);
        text2color.put("palevioletred",palevioletred);
        text2color.put("papayawhip",papayawhip);
        text2color.put("peachpuff",peachpuff);
        text2color.put("peru",peru);
        text2color.put("plum",plum);
        text2color.put("powderblue",powderblue);
        text2color.put("purple",purple);
        text2color.put("rosybrown",rosybrown);
        text2color.put("royalblue",royalblue);
        text2color.put("saddlebrown",saddlebrown);
        text2color.put("salmon",salmon);
        text2color.put("sandybrown",sandybrown);
        text2color.put("seagreen",seagreen);
        text2color.put("seashell",seashell);
        text2color.put("sienna",sienna);
        text2color.put("silver",silver);
        text2color.put("skyblue",skyblue);
        text2color.put("slateblue",slateblue);
        text2color.put("slategray",slategray);
        text2color.put("slategrey",slategrey);
        text2color.put("snow",snow);
        text2color.put("springgreen",springgreen);
        text2color.put("steelblue",steelblue);
        text2color.put("tan",tan);
        text2color.put("teal",teal);
        text2color.put("thistle",thistle);
        text2color.put("tomato",tomato);
        text2color.put("turquoise",turquoise);
        text2color.put("violet",violet);
        text2color.put("wheat",wheat);
        text2color.put("whitesmoke",whitesmoke);
        text2color.put("yellowgreen",yellowgreen);
        //
        text2color.put("black",Color.BLACK);
        text2color.put("red",Color.RED);
        text2color.put("cyan",Color.CYAN);
        text2color.put("blue",Color.BLUE);
        text2color.put("orange",Color.ORANGE);
        text2color.put("yellow",Color.YELLOW);
        text2color.put("gray",Color.GRAY);
        text2color.put("lightgray",Color.LIGHT_GRAY);
        text2color.put("darkgray",Color.DARK_GRAY);
        text2color.put("green",Color.GREEN);
        text2color.put("white",Color.WHITE);
        text2color.put("magenta",Color.MAGENTA);
        text2color.put("pink",Color.PINK);
        }
    
    public Set<String> getColorNames()
    	{
    	return Collections.unmodifiableSet(this.text2color.keySet());
    	}
        
    
    /**
     * Inverse RGB of a color
     * @param c
     * @return
     */
    public Color negative(final Color c)
    	{
    	return new Color(
    		255-c.getRed(),
    		255-c.getGreen(),
    		255-c.getBlue()
    		);
    	}
    
    /** return a gray color */
    public Color gray(float f)
    	{
    	return new Color(f,f,f);
    	}
    
    /**
     * Choose a gradient color between two colors
     * @param first the first Color
     * @param second the second Color
     * @param ratio the fraction of the first/second colors
     * @return the gradient color
     */
    public Color between(Color first,Color second,double ratio)
    	{
    	if(ratio<0 || ratio>1.0) throw new IllegalArgumentException("0<=ratio<=1 but ratio="+ratio);
    	return new Color(
    			(int)(first.getRed()+ (second.getRed()-first.getRed())*ratio),
    			(int)(first.getGreen()+ (second.getGreen()-first.getGreen())*ratio),
    			(int)(first.getBlue()+ (second.getBlue()-first.getBlue())*ratio),
    			(int)(first.getAlpha()+ (second.getAlpha()-first.getAlpha())*ratio)
    			);
    	}
    
    /**
     * convert this color as "rgb(red,green,blue)"
     * @param c the color
     * @return null if c is null or the rgb string
     */
    public static String toRGB(final Color c)
        {
        if(c==null) return null;
        return "rgb("+c.getRed()+","+c.getGreen()+","+c.getBlue()+")";
        }
    
    /**
     * parse Color as a string.
     * String can be a nominal SVG value e.g. "blue","red"... or
     * a RGB definition such as "rgb(100,200,300)" or "#ab12cc"
     * @param c the color as a string
     * @return the Color or null if it is "none" or if it cannot convert the string
     */
    public Color parse(String c)
        {
        if(c==null) return null;
        c=c.trim().toLowerCase();
        if(c.equals("none")) return null;
        else if(c.startsWith("#"))
            {
            return new Color( Integer.valueOf(c.substring(1),16).intValue());
            }
        else if((c.startsWith("rgb(") || c.startsWith("rgba(")) && c.endsWith(")"))
            {
        	final int par = c.indexOf('(');
        	final String  tokens[]= CharSplitter.COMMA.split(c.substring(par+1, c.length()-1).trim());
            if(tokens.length<3)
                    {
                    return null;
                    }
            try
                {
                final int r= java.lang.Integer.parseInt(tokens[0]);
                final int g= java.lang.Integer.parseInt(tokens[1]);
                final int b= java.lang.Integer.parseInt(tokens[2]);
                final int a=(tokens.length>3?java.lang.Integer.parseInt(tokens[3]):255);
                return new Color(r,g,b,a);
                }
            catch(final Throwable err)
                {
                return null;
                }
            }
        final Color color= text2color.get(c);
        if(color==null) throw new IllegalArgumentException("Illegal Color:"+c);
        return color;
        }
    
    @Override
    public String toString() {
    	return "ColorUtils";
    	}
    
    /** extract the YC tag of a read */
    public static class SAMRecordColorExtractor 
    	implements Function<SAMRecord,Color>
    		{
    		/** return the Color for this read. return null on error or if there is no YC tag*/
    		@Override
    		public Color apply(final SAMRecord rec) {
    			if( rec==null) return null;
    			final Object o=rec.getAttribute(YC_TAG);
    			if(o==null || !(o instanceof String)) return null;
    			final String tag =  String.class.cast(o);
    			final int comma1 = tag.indexOf(',');
    			if( comma1 <= 0) return null;
    			final int comma2 = tag.indexOf(',',comma1+1);
    			if( comma2 == -1) return null;
    			try 
    				{
    				final int r= Integer.parseInt(tag.substring(0, comma1));
    				final int g= Integer.parseInt(tag.substring(comma1+1,comma2));
    				final int b= Integer.parseInt(tag.substring(comma2+1));
    				if(r<0 || r> 255 || g<0 || g>255 || b<0 || b>255) return null;
    				return new Color(r,g,b);
    				}
    			catch(final NumberFormatException err)
    				{
    				return null;
    				}
    			}
    		}
    	

    /** converts a java.awt.Color to a YC SAMRecord attribute value */
    public static String colorToSamAttribute(final Color color) {
    	return String.valueOf(color.getRed())+","+
    			color.getGreen()+","+
    			color.getBlue()
    			;
    	}
    
    
    public static Color average(Color...colors)
    	{
    	if(colors.length==0) throw new IllegalArgumentException("empty array");
    	double r=0,g=0,b=0,a=0;
    	for(final Color c:colors)
    		{
    		r+=c.getRed();
    		g+=c.getGreen();
    		b+=c.getBlue();
    		a+=c.getAlpha();
    		}
    	return new  Color((int)(r/colors.length),(int)(g/colors.length),(int)(b/colors.length),(int)(a/colors.length));
    	}
    
	}
