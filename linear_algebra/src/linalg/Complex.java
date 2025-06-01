package linalg;

public final class Complex implements Comparable<Complex> {
	
	//TODO right now the inverse of zero is a mythical number called NaNNaNi
	
	public final double re;
	public final double im;
	
	
	public Complex(double x, double y) {
		this.re = x;
		this.im = y;
	}
	
	public Complex(double x) {
		this.re = x;
		this.im = 0;
	}
	
	public static Complex cis(double theta) {
		return new Complex(Math.cos(theta), Math.sin(theta));
	}
	
	public static Complex rcis(double r, double theta) {
		return scale(r, cis(theta));
	}
	
	
	// defaults
	
	public static final Complex ZERO = new Complex(0,0);
	public static final Complex ONE = new Complex(1,0);
	public static final Complex I = new Complex(0,1);
	public static final Complex i = I;
	public static final Complex OMEGA3 = new Complex(-0.5, Math.sqrt(3)/2);
	
	
	// something random
	
	public static Complex random() {
		return new Complex(Math.random(), Math.random());
		// returns a random complex number within the unit square
	}
	
	
	// repr
	
	public String toString() {
		String reString;
		String imString;
		
		// special case
		if(re==0 && im==0) {
			return "0";
		}
		
		// don't show the .0 at the end if you don't have to
		if(re == (int) re) {
			reString = "" + (int)re;
		} else {
			reString = "" + re;
		}
		
		// don't show the .0 at the end, and don't show the 1 before the i
		if(im == 1) {
			imString = "";
		} else if (im == -1) {
			imString = "-";
		} else if(im == (int) im) {
			imString = "" + (int)im;
		} else {
			imString = "" + im;
		}
		
		if(re == 0) {
			return imString + "i"; 
		}
		
		// don't show the +0i at the end if you don't have to
		if(im == 0) {
			return reString;
		}
		
		if(im > 0) {
			return reString + "+" + imString + "i";
		} else {
			return reString + imString + "i";
		}
	}
	
	
	public int compareTo(Complex that) {
		// complex numbers don't have a natural order, but having some kind of sorting order for java
		//  to use makes a lot of algorithms more efficient. For example, now it can be used with TreeMap.
		if(this.equals(that)) {
			return 0;
		}
		// it sorts by real and then by imaginary, so (2+2i) < (3+i) < (3+4i) < (4+2i)
		if(this.re > that.re) {
			return 1;
		}
		if(this.re < that.re) {
			return -1;
		}
		if(this.im > that.im) {
			return 1;
		}
		if(this.im < that.im) {
			return -1;
		}
		return 0;
	}
	
	
	
	protected static final double TOLERANCE = 0.000000000002;
	// equality (each of these methods has a tolerance to combat floating point errors)
	// the threshhold is 0.000000000002 = 2E-12, determined via trial and error
	public boolean equals(Complex that) {
		return sub(this, that).mag() < TOLERANCE;
	}
	
	public boolean nonzero() {
		return ! this.equals(ZERO);
	}
	public boolean isReal() {
		return Math.abs(im) < TOLERANCE;
	}
	
	
	// cleanup (to help combat floating point errors)
	public Complex cleanup() {
		double newRe = re;
		double newIm = im;
		
		if(Math.abs(re-Math.round(re)) < TOLERANCE) {
			newRe = Math.round(re);
		}
		if(Math.abs(im-Math.round(im)) < TOLERANCE) {
			newIm = Math.round(im);
		}
		
		return new Complex(newRe, newIm);
	}
	
	
	// unary operations
	
	public double mag2() {
		return re*re + im*im;
	}
	
	public double mag() {
		return Math.sqrt(mag2());
	}
	
	public double modulo() {
		return mag();
	}
	
	public double argument() {
		return Math.atan2(im, re);
	}
	
	public double arg() {
		return argument();
	}
	
	public Complex negate() {
		return new Complex(-re, -im);
	}
	
	public Complex conjugate() {
		return new Complex(re, -im);
	}
	
	public Complex conj() {
		return conjugate();
	}
	
	public Complex inverse() {
		return scale(conjugate(), 1/mag2());
	}
	
	public Complex exp() {
		//e^(x+iy) = e^x e^iy
		return rcis(Math.exp(re), im);
	}
	
	public Complex log() {
		return new Complex(Math.log(mag()), arg());
	}
	
	public Complex squared() {
		return mul(this,this);
	}
	
	public Complex cubed() {
		return mul(this, this, this);
	}
	
	
	// binary operation methods
	
	public Complex scale(double c) {
		return new Complex(c*re, c*im);
	}
	
	public Complex pow(double c) {
		return rcis( Math.pow(mag(), c), arg() * c );
	}
	
	
	// static operations
	
	public static double mag(Complex a) {
		return a.mag();
	}
	
	public static double abs(Complex a) {
		return a.mag();
	}
	
	public static Complex add(Complex a, Complex b) {
		return new Complex(a.re + b.re, a.im + b.im);
	}
	
	public static Complex negate(Complex a) {
		return a.negate();
	}
	
	public static Complex neg(Complex a) {
		return a.negate();
	}
	
	public static Complex sub(Complex a, Complex b) {
		return add(a, b.negate());
	}
	
	public static Complex mul(Complex a, Complex b) {
		return new Complex(a.re * b.re - a.im * b.im,
						   a.re * b.im + a.im * b.re);
	}
	
	public static Complex scale(Complex a, double c) {
		return a.scale(c);
	}
	
	public static Complex scale(double c, Complex a) {
		return a.scale(c);
	}
	
	public static Complex invert(Complex a) {
		return a.inverse();
	}
	
	public static Complex div(Complex a, Complex b) {
		return mul(a, b.inverse());
	}
	
	public static Complex div(Complex a, double b) {
		return scale(a, 1/b);
	}
	
	public static Complex exp(Complex a) {
		return a.exp();
	}
	
	public static Complex log(Complex a) {
		return a.log();
	}
	
	public static Complex cos(Complex z) {
		// cos(z) = (e^iz + e^-iz)/2
		return div( add( exp( mul(z,I) ), exp( mul(z,neg(I)) ) ), 2);
	}
	public static Complex sin(Complex z) {
		// sin(z) = (e^ix - e^-ix)/2i
		return div( sub( exp( mul(z,I) ), exp( mul(z,neg(I)) ) ), I.scale(2));
	}
	
	public static Complex arcsin(Complex z) {
		// asin(z) = -i log( iz + sqrt(1-z²) )
		return mul(neg(i), log(add( mul(z,i), sqrt(sub(ONE,z.squared())) )));
	}
	
	public static Complex pow(Complex a, double r) {
		return a.pow(r);
	}
	
	public static Complex pow(Complex a, Complex b) {
		// a^b = e^(ln(a) * b)
		return exp( mul( a.log(), b ) );
	}
	
	public static Complex sqrt(Complex a) {
		return a.pow(0.5);
	}
	
	public static Complex square(Complex a) {
		return mul(a,a);
	}
	
	public static Complex cbrt(Complex a) {
		return a.pow(1.0/3.0);
	}
	
	
	// varargs versions for sum and product
	
	public static Complex add(Complex... nums ) {
		Complex ret = ZERO;
		for(Complex a: nums) {
			ret = add(ret, a);
		}
		return ret;
	}
	
	public static Complex mul(Complex... nums ) {
		Complex ret = ONE;
		for(Complex a: nums) {
			ret = mul(ret, a);
		}
		return ret;
	}
	
	
	
	// useful method for taylor series
	// returns dᵏ/dzᵏ [sin] evaluated at z
	public static Complex sinDeriv(int k, Complex z) {
		return switch(k % 4) {
			case 0 -> sin(z);
			case 1 -> cos(z);
			case 2 -> neg(sin(z));
			case 3 -> neg(cos(z));
			default -> throw new IllegalArgumentException("Derivative order must be non-negative.");
		};
	}
	
	
}
