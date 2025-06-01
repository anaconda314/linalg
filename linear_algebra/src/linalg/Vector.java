package linalg;

public class Vector {
	
	private Complex[] array;
	public final int size;
	
	// constructors proper
	
	public Vector(Complex[] array) {
		this.size = array.length;
		this.array = array;
	}
	
	public Vector(double[] doubleArray) {
		this.size = doubleArray.length;
		
		this.array = new Complex[size];
		for(int i=0;i<size;i++) {
			this.array[i] = new Complex(doubleArray[i]);
		}
	}
	
	public Vector(int[] doubleArray) {
		this.size = doubleArray.length;
		
		this.array = new Complex[size];
		for(int i=0;i<size;i++) {
			this.array[i] = new Complex(doubleArray[i]);
		}
	}
	
	
	
	
	// zero vector constructor
	public static Vector ZERO(int size) {
		double[] array = new double[size];
		for(int i=0;i<size;i++) {
			array[i] = 0;
		}
		return new Vector(array);
	}
	
	// unit vector constructor
	public static Vector UNIT(int size, int index) {
		double[] array = new double[size];
		for(int i=0;i<size;i++) {
			array[i] = (i==index) ? 1 : 0;
		}
		return new Vector(array);
	}
	
	
	
	
	
	
	
	
	
	// cast to a matrix
	public Matrix asCol() {
		return Matrix.fromCols(new Vector[] {this});
	}
	public Matrix asRow() {
		return Matrix.fromRows(new Vector[] {this});
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// equality checker
	public boolean equals(Vector that) {
		if(this.size != that.size) {
			return false;
		}
		
		for(int i=0; i<size; i++) {
			if (!this.get(i).equals(that.get(i))) {
				return false;
			}
		}
		return true;
	}
	
	public boolean nonzero() {
		return ! this.equals(ZERO(size));
	}
	
	// helpful for rref.
	// start is inclusive.
	boolean nonzeroPast(int start) {
		for(int i=start; i<size; i++) {
			if(get(i).nonzero()) {
				return true;
			}
		}
		return false;
	}
	
	
	//array access
	public Complex[] asArray() {
		return array;
	}
	// getter
	public Complex get(int i) {
		return array[i];
	}
	
	
	// toString
	public String toString() {
		String ret = "[";
		for(int i=0;i<size;i++) {
			ret += array[i].toString();
			ret += (i==size-1) ? "]" : ", ";
		}
		return ret;
	}
	
	
	
	
	
	
	
	
	// vector operations
	// add, scale, negate, subtract, dot product
	
	// vector addition
	public static Vector add(Vector a, Vector b) {
		if(! (a.size == b.size) ) {
			throw new DimensionMismatchException(
					"You cannot add/subtract two vectors of different sizes.");
		}
		
		Complex[] newArray = new Complex[a.size];
		
		for(int i=0; i<a.size; i++) {
			newArray[i] = Complex.add(a.get(i), b.get(i));
		}
		return new Vector(newArray);
	}
	
	
	// scalar multiplication
	public Vector scale(Complex c) {
		Complex[] newArray = new Complex[size];
		
		for(int i=0; i<size; i++) {
			newArray[i] = Complex.mul(c, get(i));
		}
		return new Vector(newArray);
	}
	
	public Vector scale(double c) {
		return scale(new Complex(c));
	}
	public static Vector scale(Complex c, Vector v) {
		return v.scale(c);
	}
	public static Vector scale(Vector v, Complex c) {
		return v.scale(c);
	}
	public static Vector scale(double c, Vector v) {
		return v.scale(c);
	}
	public static Vector scale(Vector v, double c) {
		return v.scale(c);
	}
	
	//negation and subtraction
	public Vector negate() {
		return this.scale(-1);
	}
	public static Vector neg(Vector v) {
		return v.negate();
	}
	public static Vector sub(Vector a, Vector b) {
		return add(a,neg(b));
	}
	
	
	// and finally, dot product!
	public static Complex dot(Vector a, Vector b) {
		if(! (a.size == b.size) ) {
			throw new DimensionMismatchException(
					"You cannot dot two vectors of different sizes.");
		}
		
		Complex prod = Complex.ZERO;
		for(int i=0;i<a.size;i++) {
			prod = Complex.add(prod, Complex.mul(a.get(i),b.get(i)));
		}
		return prod;
	}
	
	
	
	
	
	
	
	
	
	
	// will come in handy at some point
	boolean isUnitVector() {
		// we need to make sure there is exactly one 1, and that there are no other nonzero numbers.
		int onesCount = 0;
		for(int i=0; i<size; i++) {
			if(get(i).equals(Complex.ONE)) {
				onesCount++;
			} else if( get(i).nonzero() ) {
				return false;
			}
		}
		
		return onesCount==1;
	}
	

}
