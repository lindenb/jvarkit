package com.github.lindenb.jvarkit.tools.treepack;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.awt.geom.Rectangle2D;

public class TreePacker
	{
	private enum Orientation { VERTICAL, HORIZONTAL};
    private enum Direction { ASCENDING, DESCENDING}
	
	private Comparator<TreePack> comparator=new Comparator<TreePack>()
				{
				public int compare(TreePack pack1, TreePack pack2)
						{
						double w1=pack1.getWeight();
						double w2=pack2.getWeight();
						if(w1<w2) return 1;
						if(w1>w2) return -1;
						return 0;
						}
				};
    
	public void layout(List<TreePack> items,final Rectangle2D bounds)
    	{
        layout(sortDescending(items),0,items.size()-1,bounds);
    	}
    
	private double sum(List<TreePack> items, int start, int end)
	    {
	    double sum=0;
	    while(start<=end)//yes <=
        	{
        	sum+=items.get(start++).getWeight();
        	}
	    return sum;
	    }
	
	private List<TreePack> sortDescending(List<TreePack> items)
	    {
	    List<TreePack> L=new ArrayList<TreePack>(items);
	    Collections.sort(L,comparator);
	    return L;
	    }
	
	private void layout(List<TreePack> items, int start, int end,
			final Rectangle2D bounds)
    {
        if (start>end) return;
            
        if (end-start<2)
        {
            layoutBest(items,start,end,bounds);
            return;
        }
        
        double x=bounds.getX(), y=bounds.getY(), w=bounds.getWidth(), h=bounds.getHeight();
        
        double total=sum(items, start, end);
        int mid=start;
        double a=items.get(start).getWeight()/total;
        double b=a;
        
        if (w<h)
        {
            // height/width
            while (mid<=end)
            {
                double aspect=normAspect(h,w,a,b);
                double q=items.get(mid).getWeight()/total;
                if (normAspect(h,w,a,b+q)>aspect) break;
                mid++;
                b+=q;
            }
            layoutBest(items,start,mid,new Rectangle2D.Double(x,y,w,h*b));
            layout(items,mid+1,end,new Rectangle2D.Double(x,y+h*b,w,h*(1-b)));
        }
        else
        {
            // width/height
            while (mid<=end)
            {
                double aspect=normAspect(w,h,a,b);
                double q=items.get(mid).getWeight()/total;
                if (normAspect(w,h,a,b+q)>aspect) break;
                mid++;
                b+=q;
            }
           layoutBest(items,start,mid,new Rectangle2D.Double(x,y,w*b,h));
           layout(items,mid+1,end,new Rectangle2D.Double(x+w*b,y,w*(1-b),h));
        }
        
    }
    
    private double aspect(double big, double small, double a, double b)
    {
        return (big*b)/(small*a/b);
    }
    
    private double normAspect(double big, double small, double a, double b)
    {
        double x=aspect(big,small,a,b);
        if (x<1) return 1/x;
        return x;
    }

    private void layoutBest(List<TreePack> items, int start, int end, final Rectangle2D bounds)
	    {
	    sliceLayout(
	    		items,start,end,bounds,
	            bounds.getWidth()>bounds.getHeight() ? Orientation.HORIZONTAL : Orientation.VERTICAL,
	            Direction.ASCENDING);
	    }
    
    

    private  void sliceLayout(List<TreePack> items, int start, int end, final Rectangle2D bounds, Orientation orientation, Direction order)
        {
            double total=sum(items, start, end);
            double a=0;
            boolean vertical=orientation==Orientation.VERTICAL;
           
            for (int i=start; i<=end; i++)
            {
            	Rectangle2D.Double r=new Rectangle2D.Double();
                double b=items.get(i).getWeight()/total;
                if (vertical)
                {
                    r.x=bounds.getX();
                    r.width=bounds.getWidth();
                    if (order==Direction.ASCENDING)
                        r.y=bounds.getY()+bounds.getHeight()*a;
                    else
                        r.y=bounds.getY()+bounds.getHeight()*(1-a-b);
                    r.height=bounds.getHeight()*b;
                }
                else
                {
                    if (order==Direction.ASCENDING)
                        r.x=bounds.getX()+bounds.getWidth()*a;
                    else
                        r.x=bounds.getX()+bounds.getWidth()*(1-a-b);
                    r.width=bounds.getWidth()*b;
                    r.y=bounds.getY();
                    r.height=bounds.getHeight();
                }
     
                items.get(i).setBounds(r);
                a+=b;
            }
        }

	}
