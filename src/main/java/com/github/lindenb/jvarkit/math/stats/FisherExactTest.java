package com.github.lindenb.jvarkit.math.stats;

import java.util.function.DoubleSupplier;

/* http://lh3lh3.users.sourceforge.net/fisher.shtml 
 * https://github.com/molgenis/systemsgenetics/blob/master/genetica-libraries/src/main/java/umcg/genetica/math/stats/FisherExactTest.java
 * 
 * 
 * */

public class FisherExactTest implements DoubleSupplier,Comparable<FisherExactTest> {
    private double left;
    private double right;
    private double twotail;
    private double sleft;
    private double sright;
    private double sless;
    private double slarg;
    private int sn11;
    private int sn1_;
    private int sn_1;
    private int sn;
    private double sprob;
    private int n11_;
    private int n12_;
    private int n21_;
    private int n22_;
    
    private FisherExactTest() {
    }
    
    public static FisherExactTest compute(final int array[]) {
	if(array==null || array.length!=4) throw new IllegalArgumentException("array null or length!=4");
    return compute(array[0],array[1],array[2],array[3]);
    }
    
    public static FisherExactTest compute(int n11, int n12, int n21, int n22) {
    	final FisherExactTest batondecolin =new FisherExactTest();
    	batondecolin._fisher(n11, n12, n21, n22);
    	return batondecolin;
    }
    
    private double _fisher(int n11, int n12, int n21, int n22) {
        n11_ = n11;
        n12_ = n12;
        n21_ = n21;
        n22_ = n22;
        if(n11_ < 0)
            n11_ *= -1;
        if(n12_ < 0)
            n12_ *= -1;
        if(n21_ < 0)
            n21_ *= -1;
        if(n22_ < 0)
            n22_ *= -1;
        int n1_ = n11_ + n12_;
        int n_1 = n11_ + n21_;
        int n = n11_ + n12_ + n21_ + n22_;
        exact(n11_, n1_, n_1, n);
        left = sless;
        right = slarg;
        twotail = sleft + sright;
        if(twotail > 1.0D)
            twotail = 1.0D;
        return twotail;
    }
    
    public double getFisherLeftTail() {
        return left;
    }
    
    public double getFisherRightTail() {
        return right;
    }
    
    public double calculateFisherTwoTail() {
        return twotail;
    }
    
    private static double lngamm(int z) {
        double x = 0.0D;
        x += 1.6594701874084621E-07D / (double)(z + 7);
        x += 9.9349371139307475E-06D / (double)(z + 6);
        x -= 0.1385710331296526D / (double)(z + 5);
        x += 12.50734324009056D / (double)(z + 4);
        x -= 176.61502914983859D / (double)(z + 3);
        x += 771.32342877576741D / (double)(z + 2);
        x -= 1259.1392167222889D / (double)(z + 1);
        x += 676.52036812188351D / (double)z;
        x += 0.99999999999951827D;
        return (Math.log(x) - 5.5810614667953278D - (double)z) + ((double)z - 0.5D) * Math.log((double)z + 6.5D);
    }
    
    private static double lnfact(int n) {
        if(n <= 1)
            return 0.0D;
        else
            return lngamm(n + 1);
    }
    
    private static double lnbico(int n, int k) {
        return lnfact(n) - lnfact(k) - lnfact(n - k);
    }
    
    private static double hyper_323(int n11, int n1_, int n_1, int n) {
        return Math.exp((lnbico(n1_, n11) + lnbico(n - n1_, n_1 - n11)) - lnbico(n, n_1));
    }
    
    private double hyper(int n11) {
        return hyper0(n11, 0, 0, 0);
    }
    
    private double hyper0(int n11i, int n1_i, int n_1i, int ni) {
        if((n1_i == 0) & (n_1i == 0) & (ni == 0)) {
            if(n11i % 10 != 0) {
                if(n11i == sn11 + 1) {
                    sprob = sprob * (((double)sn1_ - (double)sn11) / (double)n11i) * (((double)sn_1 - (double)sn11) / (((double)n11i + (double)sn) - (double)sn1_ - (double)sn_1));
                    sn11 = n11i;
                    return sprob;
                }
                if(n11i == sn11 - 1) {
                    sprob = sprob * ((double)sn11 / ((double)sn1_ - (double)n11i)) * ((((double)sn11 + (double)sn) - (double)sn1_ - (double)sn_1) / ((double)sn_1 - (double)n11i));
                    sn11 = n11i;
                    return sprob;
                }
            }
            sn11 = n11i;
        } else {
            sn11 = n11i;
            sn1_ = n1_i;
            sn_1 = n_1i;
            sn = ni;
        }
        sprob = hyper_323(sn11, sn1_, sn_1, sn);
        return sprob;
    }
    
    private double exact(int n11, int n1_, int n_1, int n) {
        int max = n1_;
        if(n_1 < max)
            max = n_1;
        int min = (n1_ + n_1) - n;
        if(min < 0)
            min = 0;
        if(min == max) {
            sless = 1.0D;
            sright = 1.0D;
            sleft = 1.0D;
            slarg = 1.0D;
            return 1.0D;
        }
        double prob = hyper0(n11, n1_, n_1, n);
        sleft = 0.0D;
        double p = hyper(min);
        int i;
        for(i = min + 1; p < 0.99999998999999995D * prob; i++) {
            sleft += p;
            p = hyper(i);
        }
        
        i--;
        if(p < 1.0000000099999999D * prob)
            sleft += p;
        else
            i--;
        sright = 0.0D;
        p = hyper(max);
        int j;
        for(j = max - 1; p < 0.99999998999999995D * prob; j--) {
            sright += p;
            p = hyper(j);
        }
        
        j++;
        if(p < 1.0000000099999999D * prob)
            sright += p;
        else
            j++;
        if(Math.abs(i - n11) < Math.abs(j - n11)) {
            sless = sleft;
            slarg = (1.0D - sleft) + prob;
        } else {
            sless = (1.0D - sright) + prob;
            slarg = sright;
        }
        return prob;
    }
    
	@Override
	public double getAsDouble() {
		return this.calculateFisherTwoTail();
		}
    @Override
    public String toString() {
    	return String.valueOf(getAsDouble());
    	}
    
    @Override
    public int compareTo(final FisherExactTest o) {
    	return Double.valueOf( this.getAsDouble()).compareTo(Double.valueOf(o.getAsDouble()));
    	}
    
    public static void main(String[] args) {
		if(args.length!=4) {
			System.err.println("Fisher: A1 B1 A2 B2");
			return;
		}
		System.out.println(
				FisherExactTest.compute(
					Integer.parseInt(args[0]),
					Integer.parseInt(args[1]),
					Integer.parseInt(args[2]),
					Integer.parseInt(args[3])
				));
	}
    
	}
