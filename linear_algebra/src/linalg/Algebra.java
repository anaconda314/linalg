package linalg;

import static linalg.Complex.*;



public final class Algebra {
	
	
	// class of static methods
	
	
	// polynomial solvers
	
	public static Complex[] solveQuadratic(Complex a, Complex b, Complex c) {
		// [b + sqrt(b²-4ac)] / 2a
		
		Complex discriminant = sub( mul(b,b), scale(4, mul(a,c)));
		
		return new Complex[] {
				div(add(neg(b), sqrt(discriminant)), scale(2,a)).cleanup(),
				div(sub(neg(b), sqrt(discriminant)), scale(2,a)).cleanup()
		};
	}
	
	
	public static Complex[] solveCubic(Complex a, Complex b, Complex c, Complex d) {
		// special thanks to https://youtu.be/zHO3YVC2T4E?si=jg9pe0rHI5M1hfFB
		//  for helping me finally understand this.
		
		// WLOG assume that a = 1
		if(! a.equals(ONE)) {
			return solveCubic(ONE, div(b, a), div(c, a), div(d,a));
		}
		
		if(b.nonzero()) {
			// depress the cubic and solve Y-Y
			
			// p = c - b²/3
			Complex p = sub(c, div(b.squared(), 3));
			// q = d - cb/3 + 2b³/27
			Complex q = add(d, div(mul(c,b),-3), scale(b.cubed(), 2.0/27.0) );
			
			// we now have a depressed cubic 0 = t³ + pt + q
			Complex[] solutions = solveDepressedCubic(p,q);
			for(int i=0; i<3; i++) {
				// x = t - b/3
				solutions[i] = sub(solutions[i], div(b,3));
			}
			return solutions;
			
		} else {
			// the cubic is already depressed T-T
			return solveDepressedCubic(c,d);
		}
	}
	
	public static Complex[] solveDepressedCubic(Complex p, Complex q) {
		
		// the solution to x³ + px + q = 0
		
		if(p.equals(ZERO)) {
			// x³ + q = 0
			// x³ = -q
			// x = cbrt(-q)
			
			Complex soln1 = cbrt(neg(q));
			
			return new Complex[] {soln1, mul(soln1,OMEGA3), mul(soln1,OMEGA3.conj())};
		}
		
		Complex Q = div(q, 2);
		Complex P = div(p, 3);
		// now the cubic is in the form x³ + 3Px + 2Q which saves us division down the line
		
		Complex delta = sqrt( add( Q.squared(), P.cubed() ) );
		
		Complex A1 = cbrt( sub(delta, Q)  );
		Complex A2 = mul( A1, OMEGA3 );
		Complex A3 = mul( A1, OMEGA3.conj() );
		
		Complex B1 = div( neg(P) , A1);
		Complex B2 = div( neg(P) , A2);
		Complex B3 = div( neg(P) , A3);
		
		return new Complex[] {
				add(A1, B1).cleanup(),
				add(A2, B2).cleanup(),
				add(A3, B3).cleanup()
		};
		
	}
	
	
	
	
	public static Complex[] solveQuartic(Complex m, Complex a, Complex b, Complex c, Complex d) {
		// the cubic solver youtuber guy had a video on solving the quartic, but it had mistakes
		//  with extraneous solutions. A LOT of them. This solution is from wikipedia:
		//  https://en.wikipedia.org/wiki/Quartic_equation#Summary_of_Ferrari's_method
		
		// WLOG assume that a = 1
		if(! m.equals(ONE)) {
			return solveQuartic(ONE, div(a, m), div(b, m), div(c,m), div(d,m));
		}
		
		if(a.nonzero()) {
			// depress the quartic and solve
			
			// p = b - (3/8)a²
			Complex p = sub(b, scale(a.squared(), 3.0/8.0));
			// q = c - (1/16)a³ - (1/2)ap
			Complex q = add(c, div(a.cubed(),-16), div(mul(a,p), -2) );
			// r = d - (1/256)a⁴ - (1/16)a²p - (1/4)aq
			Complex r = add(d, div(a.pow(4),-256), div(mul(a.squared(),p),-16), div(mul(a,q),-4) );
			
			// we now have a depressed quartic 0 = y⁴ + py² + qy + r
			//  where y = x + a/4
			Complex[] solutions = solveDepressedQuartic(p,q,r);
			for(int i=0; i<4; i++) {
				// x = y - a/4
				solutions[i] = sub(solutions[i], div(a,4));
			}
			return solutions;
			
		} else {
			// the cubic is already depressed T-T
			return solveDepressedQuartic(b,c,d);
		}
		
	}
	
	public static Complex[] solveDepressedQuartic(Complex a, Complex b, Complex c) {
		// a, b, and c are repurposed, following wikipedia's notation.
		// x⁴ + ax² + bx + c = 0
		
		if(c.equals(ZERO)) {
			// This case easily reduces to a cubic.
			// x⁴ + ax² + bx = 0
			// x (x³ + ax + b) = 0
			Complex[] solns = solveDepressedCubic(a,b);
			return new Complex[] {ZERO, solns[0], solns[1], solns[2]};
		}
		if(b.equals(ZERO)) {
			// This case easily reduces to a quadratic.
			// x⁴ + ax² + c = 0
			// (x²)² + a(x²) + c = 0
			Complex[] solns = solveQuadratic(ONE,a,c);
			return new Complex[] {sqrt(solns[0]), neg(sqrt(solns[0])), sqrt(solns[1]), neg(sqrt(solns[1]))};
		}
		
		Complex P = sub( div(a.squared(),-12), c );
		Complex Q = add( div(a.cubed(),-108), div(mul(a,c),3), div(b.squared(),-8) );
		Complex R = add( div(Q,-2), sqrt( add(div(Q.squared(),4), div(P.cubed(),27) )) );
		Complex U = cbrt(R);
		Complex y = add( a.scale(-5.0/6.0), (U.equals(ZERO)) ? neg(cbrt(Q)) : sub(U, div(P, U.scale(3))));
		Complex W = sqrt(add(a,y.scale(2)));
		
		return new Complex[] {
			add(W, sqrt(neg(add( a.scale(3), y.scale(2), div(b,W).scale(2) ))) ).scale(0.5),
			sub(W, sqrt(neg(add( a.scale(3), y.scale(2), div(b,W).scale(2) ))) ).scale(0.5),
			add(neg(W), sqrt(neg(add( a.scale(3), y.scale(2), div(b,W).scale(-2) ))) ).scale(0.5),
			sub(neg(W), sqrt(neg(add( a.scale(3), y.scale(2), div(b,W).scale(-2) ))) ).scale(0.5)
		};	
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	public static Complex[] DurandKernerFindUniqueRoots(Polynomial p) {
		Complex[] roots = new Complex[p.degree];
		for(int i=0; i<p.degree; i++) {
			roots[i] = random();
		}
		
		int iterations = 0;
		while(iterations < 1000) {
			Complex[] newRoots = roots.clone();
			for(int i=0; i<p.degree; i++) {
				
				Complex r = roots[i];
				// r_new = r - f(r)/PROD[r-s] for each other root s
				// concurrent modification is okay
				Complex prod = ONE;
				for(int j=0; j<p.degree; j++) {
					prod = mul(prod, j==i ? ONE : sub(r,newRoots[j]));
				}
				newRoots[i] = sub(r, div(p.substitute(r), prod));
			}
			
			boolean weDoneYet = true;
			for(int i=0; i<p.degree; i++) {
				if( mag(sub(newRoots[i],roots[i])) >= TOLERANCE/10 ) {
					weDoneYet = false;
				}
			}
			if(weDoneYet) break;
			
			roots = newRoots;
		
		}
		
		return roots;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// polynomial solvers (real-number versions)
	
	public static Complex[] solveQuadratic(double a, double b, double c) {
		return solveQuadratic(new Complex(a), new Complex(b), new Complex(c));
	}
	public static Complex[] solveCubic(double a, double b, double c, double d) {
		return solveCubic(new Complex(a), new Complex(b), new Complex(c), new Complex(d));
	}
	public static Complex[] solveCubic(double a, double b, double c, double d, double e) {
		return solveQuartic(new Complex(a), new Complex(b), new Complex(c), new Complex(d), new Complex(e));
	}
	
	// polynomial solver (polynomial versions)
	
	public static Complex[] solveQuadratic(Polynomial p) {
		if(p.degree != 2) {
			throw new IllegalArgumentException("The quadratic solver expected a quadratic.");
		}
		return solveQuadratic(p.term(2),p.term(1),p.term(0));
	}
	public static Complex[] solveCubic(Polynomial p) {
		if(p.degree != 3) {
			throw new IllegalArgumentException("The cubic solver expected a cubic.");
		}
		return solveCubic(p.term(3),p.term(2),p.term(1),p.term(0));
	}
	public static Complex[] solveQuartic(Polynomial p) {
		if(p.degree != 4) {
			throw new IllegalArgumentException("The quartic solver expected a quartic.");
		}
		return solveQuartic(p.term(4),p.term(3),p.term(2),p.term(1),p.term(0));
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// ultimate solver
	
	public static Complex[] solve(Polynomial p) {
		// calls the real solve method and then cleans up the output.
		Complex[] ret = solveMessy(p);
		
		for(int i=0; i<p.degree; i++) {
			Complex z = ret[i];
			
			double newRe = z.re;
			double newIm = z.im;
			if( Math.abs(z.re - Math.round(z.re)) < 0.00002) {
				newRe = Math.round(z.re);
			}
			if( Math.abs(z.im - Math.round(z.im)) < 0.00002) {
				newIm = Math.round(z.im);
			}
			
			ret[i] = new Complex(newRe, newIm);
		}
		return ret;
	}
	
	public static Complex[] solveMessy(Polynomial p) {
		if(p.degree == -1) {
			throw new IllegalArgumentException("The zero polynomial has infinite roots.");
		}
		if(p.degree == 0) {
			return new Complex[] {};
		}
		if(p.degree == 1) {
			// p(x) = ax + b
			// root is -b/a
			return new Complex[] { div(p.constantTerm(),p.linearTerm()).negate() };
		}
		if(p.degree == 2) {
			return solveQuadratic(p);
		}
		if(p.degree == 3) {
			return solveCubic(p);
		}
		if(p.degree == 4) {
			return solveQuartic(p);
		}
		
		if(p.constantTerm().equals(ZERO)) {
			Complex[] solnsMissingOne = solve(Polynomial.div(p,Polynomial.X));
			Complex[] solns = new Complex[p.degree];
			for(int i=0; i<p.degree-1; i++) {
				solns[i] = solnsMissingOne[i];
			}
			solns[p.degree-1] = ZERO;
			return solns;
		}
		
		if(p.isEven()) {
			Complex[] coeffs = new Complex[p.degree/2 + 1];
			for(int i=0; i<=p.degree; i+=2) {
				coeffs[i/2] = p.term(i);
			}
			Polynomial pInXSquared = new Polynomial(coeffs);
			
			Complex[] solns = new Complex[p.degree];
			int count = 0;
			for(Complex solnSquared: solve(pInXSquared)) {
				solns[count] = sqrt(solnSquared);
				solns[count+1] = neg(sqrt(solnSquared));
				count+=2;
			}
			return solns;
		}
		
		// if there are multiple-roots, we can find them.
		
		Polynomial gcd = Polynomial.GCD(p,p.diff());
		if(gcd.degree > 0) {
			
			Complex[] roots = new Complex[p.degree];
			Complex[] repeatRoots = solve(gcd);
			
			for(int i=0; i<gcd.degree; i++) {
				roots[i] = repeatRoots[i];
			}
			
			Complex[] uniqueRoots = solve(Polynomial.div(p,gcd));
			
			for(int i=gcd.degree; i<p.degree; i++) {
				roots[i] = uniqueRoots[i-gcd.degree];
			}
			
			return roots;
		}
		
		// now we know there are no repeat roots, so we can use Durand-Kerner
		
		return DurandKernerFindUniqueRoots(p);
		
		
		// We no longer have to throw an exception! All polynomials can now be solved!
		/*else {
			throw new UnsupportedOperationException(
					"Solutions for degree 5+ polynomials are generally impossible.");
		}*/
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
}
