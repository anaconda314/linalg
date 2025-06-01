package linalg;

public class Analysis {

	// for all the crazy functions I will need in order to take arcsin of a matrix.
	
	public static int factorial(int n) {
		if(n==0) {
			return 1;
		}
		if(n<0) {
			throw new ArithmeticException("factorial of a negative number");
		}
		return n*factorial(n-1);
	}
	
	public static double gamma(double n) {
		if(Math.floor(n*2) != n*2) {
			throw new UnsupportedOperationException(
					"The gamma function is only supported for integers and half-integers.");
		}
		if(n==1) {
			return 1;
		}
		if(n==0.5) {
			return Math.sqrt(Math.PI);
		}
		if(n<0) {
			throw new ArithmeticException("The gamma function is only supported for nonnegative numbers.");
		}
		return (n-1)*gamma(n-1);
	}
	
	
	
	public static Polynomial jacobiPolynomial(int n, int alpha, int beta) {
		if(alpha < -1 || beta < -1) {
			throw new IllegalArgumentException("alpha and beta for the Jacobi Polynomial must be >= -1.");
		}
		if(n < 0) {
			throw new IllegalArgumentException("n for the Jacobi Polynomial must be nonnegative.");
		}
		// P(alpha,beta,n)(x) = (-1)ⁿ/2ⁿn! * (1-x)⁻ᵅ(1+x)⁻ᵝ * dⁿ/dxⁿ [(1-x)ᵅ⁺ⁿ(1+x)ᵝ⁺ⁿ]
		
		Polynomial aa = (new Polynomial(-1,1));
		Polynomial bb = (new Polynomial(1,1));
		
		return
			Polynomial.divExact(
					Polynomial.mul( aa.pow(alpha+n), bb.pow(beta+n) ).diff(n),
					Polynomial.mul( aa.pow(alpha), bb.pow(beta))
			).scale( Math.pow(-0.5, n) / factorial(n) );		
	}
	
	public static Polynomial legendrePolynomialHelper(int n) {
		// https://mathworld.wolfram.com/LegendrePolynomial.html
		if(n==0) {
			return Polynomial.ONE;
		}
		if(n==1) {
			return Polynomial.X;
		}
		
		Complex c1 = Complex.div(
			Polynomial.mul(
				Polynomial.X,
				legendrePolynomialHelper(n-1).pow(2)
			).integrate(-1,1),
			legendrePolynomialHelper(n-1).pow(2).integrate(-1,1)
		);
		
		Complex c2 = Complex.div(
			legendrePolynomialHelper(n-1).pow(2).integrate(-1,1),
			legendrePolynomialHelper(n-2).pow(2).integrate(-1,1)
		);
		
		return Polynomial.sub(
			Polynomial.mul(
				legendrePolynomialHelper(n-1),
				Polynomial.sub(Polynomial.X,new Polynomial(c1))
			),
			Polynomial.scale(
				legendrePolynomialHelper(n-2),
				c2
			)
		).cleanup();
			
	}
	
	public static Polynomial legendrePolynomial(int n) {
		Polynomial p = legendrePolynomialHelper(n);
		return p.scale( p.substitute(1).inverse() );
	}
	
	
	public static Complex arcsinDeriv(int n, Complex z) {
		if(n<0) {
			throw new IllegalArgumentException("derivative order must be nonnegative");
		}
		if(n==0) {
			return Complex.arcsin(z);
		}
		
		Complex c = Complex.sub(Complex.ONE, z.squared()).pow(-0.5);
		
		return Complex.mul(
			Complex.I.pow(n-1),
			c.pow(n),
			new Complex(factorial(n-1)),
			legendrePolynomial(n-1).substitute(Complex.mul(Complex.I, z, c))
		);
	}
	
	
}
