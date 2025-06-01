package linalg;

class PolynomialMatrix {
	// A class that represents a matrix whose entries are symbolic polynomials.
	// The purpose of this class is just to act as a helper for calculating charpoly of a matrix.
	
	// The matrix will always be square, since rectangular polynomial matrices are not needed.
	// This class is very bare-bones, only containing the methods it absolutely needs to have to serve
	//  its one purpose.
	
	public final int size;
	private Polynomial[][] array;
	
	// here is the normal constructor.
	public PolynomialMatrix(Polynomial[][] array){
		// array SHOULD BE SQUARE
		this.size = array.length;
		if(array[0].length != size) {
			throw new DimensionMismatchException("Polynomials are only supported in square matrices.");
		};
		
		this.array = array;
	}
	// here is a constructor that essentially casts a normal complex Matrix to a PolynomialMatrix
	public PolynomialMatrix(Matrix M) {
		if(!M.isSquare()) {
			throw new DimensionMismatchException("Polynomials are only supported in square matrices.");
		}
		this.size = M.height;
		this.array = new Polynomial[size][size];
		
		for(int row=0; row<size; row++) {
			for(int col=0; col<size; col++) {
				array[row][col] = new Polynomial(M.get(row, col));
			}
		}
		
	}
	
	
	
	// here is the determinant method.
	public Polynomial determinant() {
		
		// base case since this definition is recursive
		if(size == 1) {
			return get(0,0);
		}
		
		Polynomial sum = Polynomial.ZERO;
		int sign = 1;
		for(int i=0; i<size; i++) {
			// sum +=   M[0][i] * det(M_0,i) * (-1)^i
			sum = Polynomial.add(sum, Polynomial.mul(get(0,i), minor(0,i).determinant()).scale(sign) );
			
			sign *= -1;
		} 
		
		return sum;
	}
	
	// that method needed get(i,j) and minor(i,j) so here they are
	public Polynomial get(int row, int col) {
		return array[row][col];
	}
	public PolynomialMatrix minor(int strikeRow, int strikeCol) {
		Polynomial[][] newArray = new Polynomial[size-1][size-1];
		
		for(int row=0; row<strikeRow; row++) {
			for(int col=0; col<strikeCol; col++) {
				newArray[row][col] = array[row][col];
			}
			for(int col=strikeCol+1; col<size; col++) {
				newArray[row][col-1] = array[row][col];
			}
		}
		for(int row=strikeRow+1; row<size; row++) {
			for(int col=0; col<strikeCol; col++) {
				newArray[row-1][col] = array[row][col];
			}
			for(int col=strikeCol+1; col<size; col++) {
				newArray[row-1][col-1] = array[row][col];
			}
		}
		return new PolynomialMatrix(newArray);
	}
	
	
	
	
	
	
	
	
	// identity and scaling will be useful for λI
	public static PolynomialMatrix IDENTITY(int size) {
		Polynomial[][] array = new Polynomial[size][size];
		for(int row=0; row<size; row++) {
			for(int col=0; col<size; col++) {
				array[row][col] = (row==col) ? Polynomial.ONE : Polynomial.ZERO;
			}
		}
		return new PolynomialMatrix(array);
	}
	public PolynomialMatrix scale(Polynomial lambda) {
		Polynomial[][] newArray = new Polynomial[size][size];
		
		for(int row=0; row<size; row++) {
			for(int col=0; col<size; col++) {
				newArray[row][col] = Polynomial.mul(lambda, get(row, col));
			}
		}
		return new PolynomialMatrix(newArray);
	}
	public PolynomialMatrix scale(double lambda) {
		return scale(new Polynomial(lambda));
	}
	
	
	
	// and finally, I'll need subtraction for A-λI
	public static PolynomialMatrix sub(PolynomialMatrix A, PolynomialMatrix B) {
		return add(A, B.scale(-1));
	}
	public static PolynomialMatrix add(PolynomialMatrix A, PolynomialMatrix B) {
		if(! (A.size == B.size) ) {
			throw new DimensionMismatchException(
					"You cannot add/subtract two matrices with different dimensions.");
		}
		
		Polynomial[][] newArray = new Polynomial[A.size][A.size];
		
		for(int row=0; row<A.size; row++) {
			for(int col=0; col<A.size; col++) {
				newArray[row][col] = Polynomial.add(A.get(row, col), B.get(row, col));
			}
		}
		return new PolynomialMatrix(newArray);
		
	}
	
	
}
