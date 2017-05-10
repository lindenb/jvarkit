/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.ParameterException;


public class ColorUtils
	{
	public static class Converter implements IStringConverter<Color> {
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
	
	
    private Map<String,Color> text2color=null;
    public ColorUtils()
		{
        text2color= new HashMap<String,Color>(150);
        text2color.put("aliceblue", new Color(240, 248, 255));
        text2color.put("antiquewhite", new Color(250, 235, 215));
        text2color.put("aquamarine", new Color(127, 255, 212));
        text2color.put("aqua", new Color( 0, 255, 255));
        text2color.put("azure", new Color(240, 255, 255));
        text2color.put("beige", new Color(245, 245, 220));
        text2color.put("bisque", new Color(255, 228, 196));
        text2color.put("black", Color.BLACK);
        text2color.put("blanchedalmond", new Color(255, 235, 205));
        text2color.put("blue", Color.BLUE);
        text2color.put("blueviolet", new Color(138, 43, 226));
        text2color.put("brown", new Color(165, 42, 42));
        text2color.put("burlywood", new Color(222, 184, 135));
        text2color.put("cadetblue", new Color( 95, 158, 160));
        text2color.put("chartreuse", new Color(127, 255, 0));
        text2color.put("chocolate", new Color(210, 105, 30));
        text2color.put("coral", new Color(255, 127, 80));
        text2color.put("cornflowerblue", new Color(100, 149, 237));
        text2color.put("cornsilk", new Color(255, 248, 220));
        text2color.put("crimson", new Color(220, 20, 60));
        text2color.put("cyan", Color.CYAN);
        text2color.put("darkblue", new Color( 0, 0, 139));
        text2color.put("darkcyan", new Color( 0, 139, 139));
        text2color.put("darkgoldenrod", new Color(184, 134, 11));
        text2color.put("darkgray", new Color(169, 169, 169));
        text2color.put("darkgreen", new Color( 0, 100, 0));
        text2color.put("darkgrey", new Color(169, 169, 169));
        text2color.put("darkkhaki", new Color(189, 183, 107));
        text2color.put("darkmagenta", new Color(139, 0, 139));
        text2color.put("darkolivegreen", new Color( 85, 107, 47));
        text2color.put("darkorange", new Color(255, 140, 0));
        text2color.put("darkorchid", new Color(153, 50, 204));
        text2color.put("darkred", new Color(139, 0, 0));
        text2color.put("darksalmon", new Color(233, 150, 122));
        text2color.put("darkseagreen", new Color(143, 188, 143));
        text2color.put("darkslateblue", new Color( 72, 61, 139));
        text2color.put("darkslategray", new Color( 47, 79, 79));
        text2color.put("darkslategrey", new Color( 47, 79, 79));
        text2color.put("darkturquoise", new Color( 0, 206, 209));
        text2color.put("darkviolet", new Color(148, 0, 211));
        text2color.put("deeppink", new Color(255, 20, 147));
        text2color.put("deepskyblue", new Color( 0, 191, 255));
        text2color.put("dimgray", new Color(105, 105, 105));
        text2color.put("dimgrey", new Color(105, 105, 105));
        text2color.put("dodgerblue", new Color( 30, 144, 255));
        text2color.put("firebrick", new Color(178, 34, 34));
        text2color.put("floralwhite", new Color(255, 250, 240));
        text2color.put("forestgreen", new Color( 34, 139, 34));
        text2color.put("fuchsia", new Color(255, 0, 255));
        text2color.put("gainsboro", new Color(220, 220, 220));
        text2color.put("ghostwhite", new Color(248, 248, 255));
        text2color.put("goldenrod", new Color(218, 165, 32));
        text2color.put("gold", new Color(255, 215, 0));
        text2color.put("gray", new Color(128, 128, 128));
        text2color.put("green", Color.GREEN);
        text2color.put("greenyellow", new Color(173, 255, 47));
        text2color.put("grey", new Color(128, 128, 128));
        text2color.put("honeydew", new Color(240, 255, 240));
        text2color.put("hotpink", new Color(255, 105, 180));
        text2color.put("indianred", new Color(205, 92, 92));
        text2color.put("indigo", new Color( 75, 0, 130));
        text2color.put("ivory", new Color(255, 255, 240));
        text2color.put("khaki", new Color(240, 230, 140));
        text2color.put("lavenderblush", new Color(255, 240, 245));
        text2color.put("lavender", new Color(230, 230, 250));
        text2color.put("lawngreen", new Color(124, 252, 0));
        text2color.put("lemonchiffon", new Color(255, 250, 205));
        text2color.put("lightblue", new Color(173, 216, 230));
        text2color.put("lightcoral", new Color(240, 128, 128));
        text2color.put("lightcyan", new Color(224, 255, 255));
        text2color.put("lightgoldenrodyellow", new Color(250, 250, 210));
        text2color.put("lightgray", new Color(211, 211, 211));
        text2color.put("lightgreen", new Color(144, 238, 144));
        text2color.put("lightgrey", new Color(211, 211, 211));
        text2color.put("lightpink", new Color(255, 182, 193));
        text2color.put("lightsalmon", new Color(255, 160, 122));
        text2color.put("lightseagreen", new Color( 32, 178, 170));
        text2color.put("lightskyblue", new Color(135, 206, 250));
        text2color.put("lightslategray", new Color(119, 136, 153));
        text2color.put("lightslategrey", new Color(119, 136, 153));
        text2color.put("lightsteelblue", new Color(176, 196, 222));
        text2color.put("lightyellow", new Color(255, 255, 224));
        text2color.put("limegreen", new Color( 50, 205, 50));
        text2color.put("lime", new Color( 0, 255, 0));
        text2color.put("linen", new Color(250, 240, 230));
        text2color.put("magenta", new Color(255, 0, 255));
        text2color.put("maroon", new Color(128, 0, 0));
        text2color.put("mediumaquamarine", new Color(102, 205, 170));
        text2color.put("mediumblue", new Color( 0, 0, 205));
        text2color.put("mediumorchid", new Color(186, 85, 211));
        text2color.put("mediumpurple", new Color(147, 112, 219));
        text2color.put("mediumseagreen", new Color( 60, 179, 113));
        text2color.put("mediumslateblue", new Color(123, 104, 238));
        text2color.put("mediumspringgreen", new Color( 0, 250, 154));
        text2color.put("mediumturquoise", new Color( 72, 209, 204));
        text2color.put("mediumvioletred", new Color(199, 21, 133));
        text2color.put("midnightblue", new Color( 25, 25, 112));
        text2color.put("mintcream", new Color(245, 255, 250));
        text2color.put("mistyrose", new Color(255, 228, 225));
        text2color.put("moccasin", new Color(255, 228, 181));
        text2color.put("navajowhite", new Color(255, 222, 173));
        text2color.put("navy", new Color( 0, 0, 128));
        text2color.put("none", null);
        text2color.put("oldlace", new Color(253, 245, 230));
        text2color.put("olivedrab", new Color(107, 142, 35));
        text2color.put("olive", new Color(128, 128, 0));
        text2color.put("orange", Color.ORANGE);
        text2color.put("orangered", new Color(255, 69, 0));
        text2color.put("orchid", new Color(218, 112, 214));
        text2color.put("palegoldenrod", new Color(238, 232, 170));
        text2color.put("palegreen", new Color(152, 251, 152));
        text2color.put("paleturquoise", new Color(175, 238, 238));
        text2color.put("palevioletred", new Color(219, 112, 147));
        text2color.put("papayawhip", new Color(255, 239, 213));
        text2color.put("peachpuff", new Color(255, 218, 185));
        text2color.put("peru", new Color(205, 133, 63));
        text2color.put("pink", Color.PINK);
        text2color.put("plum", new Color(221, 160, 221));
        text2color.put("powderblue", new Color(176, 224, 230));
        text2color.put("purple", new Color(128, 0, 128));
        text2color.put("red", Color.RED);
        text2color.put("rosybrown", new Color(188, 143, 143));
        text2color.put("royalblue", new Color( 65, 105, 225));
        text2color.put("saddlebrown", new Color(139, 69, 19));
        text2color.put("salmon", new Color(250, 128, 114));
        text2color.put("sandybrown", new Color(244, 164, 96));
        text2color.put("seagreen", new Color( 46, 139, 87));
        text2color.put("seashell", new Color(255, 245, 238));
        text2color.put("sienna", new Color(160, 82, 45));
        text2color.put("silver", new Color(192, 192, 192));
        text2color.put("skyblue", new Color(135, 206, 235));
        text2color.put("slateblue", new Color(106, 90, 205));
        text2color.put("slategray", new Color(112, 128, 144));
        text2color.put("slategrey", new Color(112, 128, 144));
        text2color.put("snow", new Color(255, 250, 250));
        text2color.put("springgreen", new Color( 0, 255, 127));
        text2color.put("steelblue", new Color( 70, 130, 180));
        text2color.put("tan", new Color(210, 180, 140));
        text2color.put("teal", new Color( 0, 128, 128));
        text2color.put("thistle", new Color(216, 191, 216));
        text2color.put("tomato", new Color(255, 99, 71));
        text2color.put("turquoise", new Color( 64, 224, 208));
        text2color.put("violet", new Color(238, 130, 238));
        text2color.put("wheat", new Color(245, 222, 179));
        text2color.put("white", Color.WHITE);
        text2color.put("whitesmoke", new Color(245, 245, 245));
        text2color.put("yellowgreen", new Color(154, 205, 50));
        text2color.put("yellow", Color.YELLOW);            
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
    public static String toRGB(Color c)
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
        else if(c.startsWith("rgb(") && c.endsWith(")"))
            {
            int index=c.indexOf(',');
            if(index==-1)
                    {
                    return null;
                    }
            try
                {
                int r= java.lang.Integer.parseInt(c.substring(4,index));
                c=c.substring(index+1);
                index=c.indexOf(',');
                if(index==-1)
                    {
                    return null;
                    }
                int g= java.lang.Integer.parseInt(c.substring(0,index));
                c=c.substring(index+1);
                index=c.indexOf(')');
                if(index==-1)
                    {
                    return null;
                    }
                int b= java.lang.Integer.parseInt(c.substring(0,index));
                if( c.substring(index+1).trim().length()!=0)
                        {
                        return null;
                        }
                return new Color(r,g,b);
                }
            catch(Exception err)
                {
                return null;
                }
            }
        final Color color= (Color)text2color.get(c);
        if(color==null) throw new IllegalArgumentException("Illegal Color:"+c);
        return color;
        }
    
    @Override
    public String toString() {
    	return "ColorUtils";
    	}
    
	}
