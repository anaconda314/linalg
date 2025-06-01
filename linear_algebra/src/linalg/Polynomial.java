package linalg;

public class Polynomial {

	protected int degree;
	private Complex[] coefficients;
	
	// constructors
	
	public Polynomial(Complex[] coefficients) {
		// coefficients go from x^0 to x^n
		coefficients = trimZeroes(coefficients);
		
		this.degree = coefficients.length-1;
		this.coefficients = coefficients;
	}
	
	public Polynomial(double... doubleCoefs) {
		// this is the more user-friendly version, mainly for use in testing, where doubleCoefs
		//  go from x^n to x^0
		
		this.degree = doubleCoefs.length-1;
		this.coefficients = new Complex[degree+1];
		for(int n=0; n<=degree; n++) {
			this.coefficients[degree-n] = new Complex(doubleCoefs[n]);
		}
		this.coefficients = trimZeroes(this.coefficients);
		this.degree = this.coefficients.length-1;
	}
	
	public Polynomial(int... doubleCoefs) {
		// just like the above but with ints
		this.degree = doubleCoefs.length-1;
		this.coefficients = new Complex[degree+1];
		for(int n=0; n<=degree; n++) {
			this.coefficients[degree-n] = new Complex(doubleCoefs[n]);
		}
		this.coefficients = trimZeroes(this.coefficients);
		this.degree = this.coefficients.length-1;
	}
	
	// java won't let me make a constructor Polynomial(Complex... coefficients), but I will at least
	//  make a constructor with a single Complex input for constants.
	public Polynomial(Complex constant) {
		if(constant.nonzero()) {
			this.degree = 0;
			this.coefficients = new Complex[] {constant};
		} else {
			this.degree = -1;
			this.coefficients = new Complex[] {};
		}
	}
	
	
	public static Polynomial fromRoots(Complex[] roots) {
		Polynomial p = ONE;
		for(Complex root: roots) {
			p = mul(p, new Polynomial( new Complex[] {root.negate(), Complex.ONE} ));
		}
		return p;
	}
	
	public static Polynomial fromRoots(double... roots) {
		Polynomial p = ONE;
		for(double root: roots) {
			p = mul(p, new Polynomial( 1, -root ));
		}
		return p;
	}
	
	
	// some default values
	
	public static final Polynomial ZERO = new Polynomial(0);
	public static final Polynomial ONE = new Polynomial(1);
	public static final Polynomial X = new Polynomial(1,0);
	
	
	
	
	
	
	// toString helper
	public static String termToString(Complex coef, int power) {
		if(coef.equals(Complex.ZERO)) {
			return "";
		}
		
		String ret = "";
		
		if(! coef.isReal()) {
			ret += "+ (" + coef + ")";
		} else if (coef.re == 1) {
			ret += "+ ";
		} else if (coef.re == -1) {
			ret += "- ";
		} else if (coef.re < 0) {
			ret += "- " + (coef.negate());
		} else {
			ret += "+ " + coef;
		}
		
		if(power == 0) {
			ret += "";
		} else if (power == 1) {
			ret += "x";
		} else {
			ret += "x^" + power;
		}
		
		if(ret.length() == 2) {
			// so that -1 and +1 appear again in the x^0 case
			ret += "1";
		}
		
		return ret;
	}
	public String toString() {
		String ret = "";
		for(int n=degree; n>=0; n--) {
			ret += termToString(coefficients[n], n) + " ";
		}
		
		if(ret.equals("")) {
			return "0";
		}
		if(ret.charAt(0) == '+') {
			return ret.substring(1);
		}
		// implied else
		return ret;
		
	}
	
	
	public Polynomial cleanup() {
		Complex[] newArray = new Complex[degree+1];
		
		for(int n=0; n<=degree; n++) {
			newArray[n] = coefficients[n].cleanup();
		}
		return new Polynomial(newArray);
		
	}
	
	public Polynomial monic() {
		if(!nonzero()) return this;
		
		Complex[] newArray = new Complex[degree+1];
		
		newArray[degree] = Complex.ONE;
		
		for(int n=degree-1; n>=0; n--) {
			newArray[n] = Complex.div(coefficients[n],coefficients[degree]);
		}
		return new Polynomial(newArray);
	}
	
	
	
	
	// getter methods
	public int degree() {
		return degree;
	}
	public Complex[] asArray() {
		return coefficients;
	}
	public Complex coefficient(int n) {
		if(n > degree || n < 0) {
			return Complex.ZERO;
		}
		return coefficients[n];
	}
	public Complex term(int n) {
		// alias for coefficient
		return coefficient(n);
	}
	public Complex constantTerm() {
		return coefficients[0];
	}
	public Complex linearTerm() {
		return coefficients[1];
	}
	public Complex quadraticTerm() {
		return coefficients[2];
	}
	public Complex cubicTerm() {
		return coefficients[3];
	}
	public Complex leadingCoefficient() {
		return coefficients[degree];
	}
	
	
	
	
	public boolean equals(Polynomial that) {
		for(int i=0; i<degree; i++) {
			if (!this.term(i).equals(that.term(i))) {
				return false;
			}
		}
		return true;
	}
	
	public boolean nonzero() {
		return this.degree != -1;
	}
	
	
	
	
	
	public boolean isOdd() {
		for(int n=0; n<=degree; n+=2) {
			if(term(n).nonzero()) {
				return false;
			}
		}
		return true;
	}
	
	public boolean isEven() {
		for(int n=1; n<=degree; n+=2) {
			if(term(n).nonzero()) {
				return false;
			}
		}
		return true;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	// operations.
	// add, scale, negate, subtract, multiply, differentiate
	
	// addition is more complicated than with vectors because you can add polynomials of different degrees.
	public static Polynomial add(Polynomial p, Polynomial q) {
		
		
		int degree = Math.max(p.degree, q.degree);
		
		Complex[] newArray = new Complex[degree+1];
		
		for(int n=0; n<=degree; n++) {
			newArray[n] = Complex.add(p.term(n), q.term(n));
		}
		
		return new Polynomial(newArray);
	}
	
	
	// scalar multiplication is the same as with vectors
	public Polynomial scale(Complex c) {
		Complex[] newArray = new Complex[degree+1];
		
		for(int n=0; n<=degree; n++) {
			newArray[n] = Complex.mul(c, term(n));
		}
		return new Polynomial(newArray);
	}
	
	public Polynomial scale(double c) {
		return scale(new Complex(c));
	}
	public static Polynomial scale(Complex c, Polynomial p) {
		return p.scale(c);
	}
	public static Polynomial scale(Polynomial p, Complex c) {
		return p.scale(c);
	}
	public static Polynomial scale(double c, Polynomial p) {
		return p.scale(c);
	}
	public static Polynomial scale(Polynomial p, double c) {
		return p.scale(c);
	}
	
	// negation and subtraction
	public Polynomial negate() {
		return this.scale(-1);
	}
	public static Polynomial neg(Polynomial p) {
		return p.negate();
	}
	public static Polynomial sub(Polynomial a, Polynomial b) {
		return add(a,neg(b));
	}
	
	
	// polynomial multiplication is a whole new thing

	public static Polynomial mul(Polynomial p, Polynomial q) {
		int degree = p.degree + q.degree;
		if(degree == -2) degree = -1;
		Complex[] newArray = new Complex[degree + 1];
		
		for(int n=0; n<=degree; n++) {
			newArray[n] = Complex.ZERO;
			
			for(int i=0; i<=n; i++) {
				newArray[n] = Complex.add(newArray[n], Complex.mul( p.term(i), q.term(n-i) ));
			}
		}
		return new Polynomial(newArray);
		
	}
	
	public Polynomial pow(int c) {
		if(c<0) {
			throw new IllegalArgumentException("cannot raise a polynomial to a negative power.");
		}
		Polynomial ret = Polynomial.ONE;
		for(int i=0; i<c; i++) {
			ret = mul(ret, this);
		}
		return ret;
	}
	
	public static Polynomial pow(Polynomial p, int c) {
		return p.pow(c);
	}
	
	
	
	// differentiation
	public Polynomial diff() {
		if(degree == -1) {
			// derivative of zero is zero.
			return this;
		}
		
		Complex[] newArray = new Complex[degree];
		//^ the size of the array is degree instead of degree+1 because differentiating reduces the degree.
		
		for(int n=1; n <= degree; n++) {
			newArray[n-1] = Complex.scale(coefficients[n], n);
		}
		
		return new Polynomial(newArray);
	}
	
	public Polynomial diff(int n) {
		// nth derivative
		if(n<0) {
			throw new IllegalArgumentException("Derivative order must be non-negative.");
		}
		
		Polynomial ret = this;
		for(int i=0; i<n; i++) {
			ret = ret.diff();
		}
		return ret;
	}
	
	
	
	
	
	
	// definite integration
	
	public Complex integrate(Complex a, Complex b) {
		// ʃ_a^b  p(x)dx
		
		Complex[] newArray = new Complex[degree+2];
		newArray[0] = Complex.ZERO;
		
		for(int n=0; n <= degree; n++) {
			newArray[n+1] = Complex.scale(coefficients[n], 1.0/n);
		}
		Polynomial antiderivative = new Polynomial(newArray);
		return Complex.sub(antiderivative.substitute(b), antiderivative.substitute(a));
	}
	public Complex integrate(double a, double b) {
		// ʃ_a^b  p(x)dx
		
		Complex[] newArray = new Complex[degree+2];
		newArray[0] = Complex.ZERO;
		
		for(int n=0; n <= degree; n++) {
			newArray[n+1] = Complex.scale(coefficients[n], 1.0/(n+1));
		}
		Polynomial antiderivative = new Polynomial(newArray);
		return Complex.sub(antiderivative.substitute(b), antiderivative.substitute(a));
	}
	
	
	
	
	
	
	
	
	// polynomial division
	public static Polynomial div(Polynomial dividend, Polynomial divisor) {
		Polynomial quotient = Polynomial.ZERO;
		while(dividend.degree >= divisor.degree) {
			Complex nextCoef = Complex.div(dividend.leadingCoefficient(), divisor.leadingCoefficient());
			Polynomial nextTerm = scale( nextCoef, X.pow(dividend.degree - divisor.degree) );
			
			quotient = add(quotient, nextTerm);
			dividend = sub(dividend, mul(divisor, nextTerm));
		}
		return quotient;
	}
	
	public static Polynomial rem(Polynomial dividend, Polynomial divisor) {
		Polynomial quotient = Polynomial.ZERO;
		while(dividend.degree >= divisor.degree) {
			Complex nextCoef = Complex.div(dividend.leadingCoefficient(), divisor.leadingCoefficient());
			Polynomial nextTerm = scale( nextCoef, X.pow(dividend.degree - divisor.degree) );
			
			quotient = add(quotient, nextTerm);
			dividend = sub(dividend, mul(divisor, nextTerm));
		}
		return dividend;
	}
	
	// polynomial division, but it throws an error if the polynomials are not divisible.
	public static Polynomial divExact(Polynomial dividend, Polynomial divisor) {
		Polynomial remainder = dividend;
		Polynomial quotient = Polynomial.ZERO;
		while(remainder.degree >= divisor.degree) {
			Complex nextCoef = Complex.div(remainder.leadingCoefficient(), divisor.leadingCoefficient());
			Polynomial nextTerm = scale( nextCoef, X.pow(remainder.degree - divisor.degree) );
			
			quotient = add(quotient, nextTerm);
			remainder = sub(remainder, mul(divisor, nextTerm));
		}
		if(remainder.equals(ZERO)) {
			return quotient;
		}else {
			System.out.println(dividend);
			System.out.println(divisor);
			System.out.println(remainder);
			throw new ArithmeticException("The polynomials do not divide.");
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	
	public Complex substitute(Complex z) {
		Complex ret = Complex.ZERO;
		for(int n=0; n<=degree; n++) {
			ret = Complex.add(ret, Complex.mul(term(n), z.pow(n)) );
		}
		return ret.cleanup();
	}
	
	public Complex substitute(double x) {
		Complex z = new Complex(x);
		
		Complex ret = Complex.ZERO;
		for(int n=0; n<=degree; n++) {
			ret = Complex.add(ret, Complex.mul(term(n), z.pow(n)) );
		}
		return ret;
	}
	
	public Complex substitute(int x) {
		Complex z = new Complex(x);
		
		Complex ret = Complex.ZERO;
		for(int n=0; n<=degree; n++) {
			ret = Complex.add(ret, Complex.mul(term(n), z.pow(n)) );
		}
		return ret;
	}
	
	
	// matrix substitution
	
	public Matrix substitute(Matrix M) {
		if(!M.isSquare()) {
			throw new DimensionMismatchException("Only square matrices can be raised to powers.");
		}
		
		Matrix ret = Matrix.ZERO(M.height);
		for(int n=0; n<=degree; n++) {
			ret = Matrix.add(ret, Matrix.scale(term(n), M.pow(n)) );
		}
		return ret;
	}
	
	// composition
	
	public Polynomial substitute(Polynomial p) {
		Polynomial ret = ZERO;
		for(int n=0; n<=degree; n++) {
			ret = Polynomial.add(ret, Polynomial.scale(term(n), p.pow(n)) );
		}
		return ret;
	}
	public static Polynomial compose(Polynomial p, Polynomial q) {
		return p.substitute(q);
	}
	
	// alias for substitute
	public Complex apply(Complex z) {return substitute(z);}
	public Complex apply(double x) {return substitute(x);}
	public Complex apply(int x) {return substitute(x);}
	public Matrix apply(Matrix M) {return substitute(M);}
	public Polynomial compose(Polynomial p) {return substitute(p);}
	
	
	
	
	
	
	
	
	
	// root checking
	
	public static boolean isRoot(Complex num, Polynomial p) {
		return p.substitute(num).cleanup().equals(Complex.ZERO);
	}
	
	
	
	
	
	
	
	
	public static Polynomial GCD(Polynomial p, Polynomial q) {
		if(p.degree < q.degree) {
			return GCD(q,p);
		}
		// WLOG, p.degree >= q.degree
		// base case
		if(! q.nonzero()) {
			return p.monic();
		}
		// the GCD is the same if we replace p by rem(p,q)
		return GCD(rem(p,q), q);
	}
	
	
	
	
	
	
	
	
	
	
	// number theory concepts that are necessary for solving a quintic
	
	public static Matrix sylvester(Polynomial p, Polynomial q) {
		Complex[][] array = new Complex[p.degree+q.degree][p.degree+q.degree];
		for(int row=0; row<q.degree; row++) {
			for(int col=0; col<array.length; col++) {
				array[row][col] = p.coefficient(p.degree+row-col);
			}
		}
		for(int row=0; row<p.degree; row++) {
			for(int col=0; col<array.length; col++) {
				array[row+q.degree][col] = q.coefficient(q.degree+row-col);
			}
		}
		return new Matrix(array);
	}
	
	public static Complex resultant(Polynomial p, Polynomial q) {
		return sylvester(p,q).determinant();
	}
	
	public Complex discriminant() {
		return Complex.div( resultant(this, this.diff()), this.leadingCoefficient().inverse());
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// useful
	// trims trailing zeroes from the list (useful for polynomials)
	private static Complex[] trimZeroes(Complex[] array) {
		int degree = array.length - 1;
		while(true) {
			if(degree == -1) {
				break;
			}else if(array[degree].nonzero()) {
				break;
			} else {
				degree--;
			}
		}
		
		Complex[] newArray = new Complex[degree+1];
		for(int n=0; n<=degree; n++) {
			newArray[n] = array[n];
		}
		return newArray;
	}
	
	
}
