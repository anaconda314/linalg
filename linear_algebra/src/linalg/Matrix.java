package linalg;

import static linalg.Analysis.factorial;

import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;

public class Matrix {
	public final int height;
	public final int width;
	
	private Complex[][] array;
	
	// proper constructors (a 2d array of either Complex, double, or int must be passed as input)
	
	public Matrix(Complex[][] array){
		// array SHOULD BE RECTANGULAR
		this.height = array.length;
		
		if(this.height == 0) this.width = 0;
		else this.width = array[0].length;
		
		this.array = array;
	}
	
	public Matrix(double[][] doubleArray) {
		this.height = doubleArray.length;
		
		if(this.height == 0) this.width = 0;
		else this.width = doubleArray[0].length;
		
		this.array = new Complex[height][width];
		for(int row=0; row<height; row++) {
			for(int col=0; col<width; col++) {
				this.array[row][col] = new Complex(doubleArray[row][col]);
			}
		}
	}
	
	public Matrix(int[][] doubleArray) {
		this.height = doubleArray.length;
		
		if(this.height == 0) this.width = 0;
		else this.width = doubleArray[0].length;
		
		this.array = new Complex[height][width];
		for(int row=0; row<height; row++) {
			for(int col=0; col<width; col++) {
				this.array[row][col] = new Complex(doubleArray[row][col]);
			}
		}
	}
	
	
	// constructor-ish methods
	
	public static Matrix fromRows(Vector[] rows) {
		int height = rows.length;
		int width = rows[0].size;
		
		Complex[][] array = new Complex[height][];
		for(int row=0; row<height; row++) {
			if (rows[row].size != width) {
				throw new DimensionMismatchException("All rows of a matrix must be the same size.");
			}
			array[row] = rows[row].asArray();
		}
		
		return new Matrix(array);
	}
	
	public static Matrix fromCols(Vector[] cols) {
		int height = cols[0].size;
		int width = cols.length;
		
		Complex[][] array = new Complex[height][width];
		for(int col=0; col<width; col++) {
			if (cols[col].size != height) {
				throw new DimensionMismatchException("All columns of a matrix must be the same size.");
			}
			for(int row=0; row<height; row++) {
				array[row][col] = cols[col].get(row);
			}
		}
		
		return new Matrix(array);
	}
	
	
	// special matrices
	
	public static Matrix IDENTITY(int size) {
		int[][] array = new int[size][size];
		for(int row=0; row<size; row++) {
			for(int col=0; col<size; col++) {
				array[row][col] = (row==col) ? 1 : 0;
			}
		}
		return new Matrix(array);
	}
	
	public static Matrix ZERO(int height, int width) {
		int[][] array = new int[height][width];
		for(int row=0; row<height; row++) {
			for(int col=0; col<width; col++) {
				array[row][col] = 0;
			}
		}
		return new Matrix(array);
	}
	
	public static Matrix ZERO(int size) {
		int[][] array = new int[size][size];
		for(int row=0; row<size; row++) {
			for(int col=0; col<size; col++) {
				array[row][col] = 0;
			}
		}
		return new Matrix(array);
	}
	
	public static Matrix diag(Complex[] evals) {
		int size = evals.length;
		Complex[][] array = new Complex[size][size];
		for(int row=0; row<size; row++) {
			for(int col=0; col<size; col++) {
				array[row][col] = (row==col) ? evals[row] : Complex.ZERO;
			}
		}
		return new Matrix(array);
	}
	
	public static Matrix jordanBlock(int size, Complex eigenvalue) {
		Complex[][] array = new Complex[size][size];
		for(int row=0; row<size; row++) {
			for(int col=0; col<size; col++) {
				array[row][col] = (col==row) ? eigenvalue : (col==row+1) ? Complex.ONE : Complex.ZERO;
			}
		}
		return new Matrix(array);
	}
	
	public static Matrix diag(Matrix[] blocks) {
		int height = 0;
		int width = 0;
		for(Matrix M: blocks) {
			height += M.height;
			width += M.width;
		}
		
		Complex[][] array = new Complex[height][width];
		int startRow = 0;
		int startCol = 0;
		for(Matrix M: blocks) {
			for(int row=0; row<M.height; row++) {
				for(int col=0; col<M.width; col++) {
					array[startRow + row][startCol + col] = M.array[row][col];
				}
			}
			startRow += M.height;
			startCol += M.width;
		}
		
		for(int row=0; row<height; row++) {
			for(int col=0; col<width; col++) {
				if(array[row][col] == null) {
					array[row][col] = Complex.ZERO;
				}
			}
		}
		
		return new Matrix(array);
		
	}
	
	
	
	
	
	
	
	
	
	
	// random matrix of integers between -5 and 5.
	// if it is square, it is likely invertible.
	public static Matrix random(int height, int width) {
		int[][] array = new int[height][width];
		
		for(int row=0; row<height; row++) {
			for(int col=0; col<width; col++) {
				array[row][col] = (int) Math.floor(11*Math.random()-5);
			}
		}
		
		return new Matrix(array);
	}
	
	// Random non-invertible (square) matrix.
	public static Matrix randomNoninvertible(int size) {
		Matrix I = IDENTITY(size);
		while(true) {
			Matrix R = random(size, size);
			if(!rref(R).equals(I)) return R;
		}
	}
	
	
	
	
	
	
	
	
	
	
	// array access
	public Complex[][] asArray() {
		return array;
	}
	
	// getter function
	public Complex get(int row, int col) {
		return array[row][col];
	}
	
	// get rows or columns of the matrix
	public Vector getRow(int row) {
		return new Vector(array[row]);
	}
	public Vector getCol(int col) {
		Complex[] newArray = new Complex[height];
		for(int row=0; row<height; row++) {
			newArray[row] = get(row,col);
		}
		return new Vector(newArray);
	}
	
	// get all rows or columns
	public Vector[] rows() {
		Vector[] ret = new Vector[height];
		for(int row=0; row<height; row++) {
			ret[row] = new Vector(array[row]);
		}
		return ret;
	}
	public Vector[] cols() {
		return transpose().rows();
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// the toString method will make the matrix look pretty, like so:
	/*[4+2i  5     3i  7
	 * 3     3+2i  7   2i
	 * 0     0     0   0 ]
	 */
	public String toString() {
		
		if(height == 0 || width == 0) {
			// o silly.
			return "[]";
		}
		
		// keep track of what each string will be and how much space each column will need
		String[][] stringArray = new String[height][width];
		int[] colLengths = new int[width];
		
		for(int row=0; row<height; row++) {
			for(int col=0; col<width; col++) {
				String theString = array[row][col].toString();
				colLengths[col] = Math.max(colLengths[col], theString.length());
				stringArray[row][col] = theString;
			}
		}
		
		String ret = "";
		
		for(int row=0; row<height; row++) {
			ret += (row==0) ? "[" : " ";
			for(int col=0; col<width; col++) {
				ret += (col==0) ? "" : "  ";
				ret += padRight(stringArray[row][col], colLengths[col]);
			}
			ret += (row==height-1) ? "]" : "\n";
		}
		
		return ret;
		
		
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	// equality checker
	
	public boolean equals(Matrix that) {
		if(this.height != that.height || this.width != that.width) {
			return false;
		}
		
		for(int row=0; row<height; row++) {
			for(int col=0; col<width; col++) {
				if (!this.get(row, col).equals(that.get(row, col))) {
					return false;
				}
			}
		}
		return true;
	}
	
	// other handy checker
	
	public boolean isSquare() {
		return height==width;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// this method is needed to clean up floating point errors, and relies on the tolerance in the
	//  equals method of Complex
	public Matrix cleanup() {
		Complex[][] newArray = array.clone();
		
		for(int row=0; row<height; row++) {
			for(int col=0; col<width; col++) {
				newArray[row][col] = array[row][col].cleanup();
			}
		}
		return new Matrix(newArray);
		
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// transpose
	
	public Matrix transpose() {
		Complex[][] newArray = new Complex[width][height];
		
		for(int row=0; row<height; row++) {
			for(int col=0; col<width; col++) {
				newArray[col][row] = array[row][col];
			}
		}
		
		return new Matrix(newArray);
	}
	
	public Matrix T() {
		return transpose();
	}
	
	
	
	// trace
	
	public Complex trace() {
		Complex ret = Complex.ZERO;
		
		if(! isSquare()) {
			throw new DimensionMismatchException("Only square matrices have traces.");
		}
		
		for(int i=0; i<height; i++) {
			ret = Complex.add(ret, get(i,i));
		}
		return ret;
	}
	
	
	
	
	// matrix addition
	
	public static Matrix add(Matrix A, Matrix B) {
		if(! (A.height == B.height && A.width == B.width) ) {
			throw new DimensionMismatchException(
					"You cannot add/subtract two matrices with different dimensions.");
		}
		
		Complex[][] newArray = new Complex[A.height][A.width];
		
		for(int row=0; row<A.height; row++) {
			for(int col=0; col<A.width; col++) {
				newArray[row][col] = Complex.add(A.get(row, col), B.get(row, col));
			}
		}
		return new Matrix(newArray);
		
	}
	
	
	
	
	// scalar multiplication
	
	public Matrix scale(Complex c) {
		Complex[][] newArray = new Complex[height][width];
		
		for(int row=0; row<height; row++) {
			for(int col=0; col<width; col++) {
				newArray[row][col] = Complex.mul(c, get(row, col));
			}
		}
		return new Matrix(newArray);
	}
	
	public Matrix scale(double c) {
		return scale(new Complex(c));
	}
	
	public static Matrix scale(Matrix M, Complex c) {
		return M.scale(c);
	}
	
	public static Matrix scale(Complex c, Matrix M) {
		return M.scale(c);
	}
	
	public static Matrix scale(Matrix M, double c) {
		return M.scale(c);
	}
	
	public static Matrix scale(double c, Matrix M) {
		return M.scale(c);
	}
	
	
	// negation and subtraction
	
	public Matrix negate() {
		return this.scale(-1);
	}
	
	public static Matrix neg(Matrix M) {
		return M.negate();
	}
	
	public static Matrix sub(Matrix A, Matrix B) {
		return add(A, neg(B));
	}
	
	
	
	
	
	
	// now for the big one: matrix multiplication
	public static Matrix mul(Matrix A, Matrix B) {
		if(A.width != B.height) {
			throw new DimensionMismatchException(
					"To multiply AB, the width of A must match the height of B.");
		}
		
		Complex[][] newArray = new Complex[A.height][B.width];
		for(int row=0; row<A.height; row++) {
			for(int col=0; col<B.width; col++) {
				newArray[row][col] = Vector.dot(A.getRow(row), B.getCol(col));
			}
		}
		return new Matrix(newArray).cleanup();
	}
	
	public static Matrix mul(Matrix... Ms ) {
		Matrix ret = IDENTITY(Ms[0].height);
		for(Matrix M: Ms) {
			ret = mul(ret, M);
		}
		return ret;
	}
	
	// multiplication of a matrix and a vector, note that this is not static
	public Vector mul(Vector v) {
		// convert vector to 1×n matrix, multiply, and convert 1×n matrix product back to a vector
		return Matrix.mul(this, v.asCol()).getCol(0);
	}
	
	public Matrix pow(int c) {
		if(!isSquare()) {
			throw new DimensionMismatchException("Only square matrices can be raised to powers.");
		}
		
		Matrix ret = IDENTITY(height);
		for(int i=0; i<c; i++) {
			ret = mul(ret,this);
		}
		return ret;
	}
	
	// if the input of .pow is an int, it computes it iteratively. If the input is double, it uses
	//  diagonalization to find the value.
	public Matrix pow(double c) {
		if(!isSquare()) {
			throw new DimensionMismatchException("Only square matrices can be raised to powers.");
		}
		
		Matrix[] d = diagonalize();
		Matrix P = d[0];
		Matrix D = d[1];
		Matrix PInverse = d[2];
		
		Complex[] lambdasPow = new Complex[height];
		for(int i=0; i<height; i++) {
			lambdasPow[i] = D.get(i,i).pow(c);
		}
		Matrix DPow = diag(lambdasPow);
		
		return mul(P, DPow, PInverse);
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// row operations
	
	// row addition, row scaling, row swapping
	
	public Matrix rowAdd(int target, int addend) {
		Vector[] rows = this.rows();
		rows[target] = Vector.add(rows[target], rows[addend]);
		return Matrix.fromRows(rows);
	}
	
	public Matrix rowAdd(int target, int addend, Complex scaleFactor) {
		Vector[] rows = this.rows();
		rows[target] = Vector.add(rows[target], rows[addend].scale(scaleFactor));
		return Matrix.fromRows(rows);
	}
	
	public Matrix rowAdd(int target, int addend, double scaleFactor) {
		Vector[] rows = this.rows();
		rows[target] = Vector.add(rows[target], rows[addend].scale(scaleFactor));
		return Matrix.fromRows(rows);
	}
	
	public Matrix rowScale(int target, Complex scaleFactor) {
		Vector[] rows = this.rows();
		rows[target] = rows[target].scale(scaleFactor);
		return Matrix.fromRows(rows);
	}
	
	public Matrix rowScale(int target, double scaleFactor) {
		Vector[] rows = this.rows();
		rows[target] = rows[target].scale(scaleFactor);
		return Matrix.fromRows(rows);
	}
	
	public Matrix rowSwap(int row1, int row2) {
		Vector[] rows = this.rows();
		Vector temp = rows[row1];
		rows[row1] = rows[row2];
		rows[row2] = temp;
		return Matrix.fromRows(rows);
	}
	
	
	
	
	
	
	
	
	
	
	
	// now for the rref helpers
	
	public Complex pivotEntry(int row) {
		for(int i=0; i<width; i++) {
			if( get(row,i).nonzero() ) {
				return get(row,i); 
			}
		}
		return Complex.ZERO;
	}
	
	public int pivotCol(int row) {
		for(int i=0; i<width; i++) {
			if( get(row,i).nonzero() ) {
				return i; 
			}
		}
		return width;
	}
	
	// this uses a not-very-strict definition of ref, where the pivot entries don't need to be 1.
	public boolean isRef() {
		for(int row=1; row<height; row++) {
			int oldPivotCol = pivotCol(row-1);
			int newPivotCol = pivotCol(row);
			
			if(newPivotCol <= oldPivotCol && newPivotCol != width) {
				// the only way two rows are allowed to have the same pivotCol is if they're entirely zeroes
				//  (which is what pivotCol==width indicates).
				// we don't have to worry about the case where new < old && new == width
				//  because that cannot happen.
				return false;
			}
		}
		return true;
	}
	
	public boolean isRref() {
		
		// keep track of where all the rows start
		int[] pivotCols = new int[height];
		for(int row=0; row<height; row++) {
			pivotCols[row] = pivotCol(row);
		}
		
		// this for loop checks that the matrix is staircase shaped
		for(int row=1; row<height; row++) {
			if(pivotCols[row] <= pivotCols[row-1] && pivotCols[row] != width) {
				// same logic as in the isRef() function
				return false;
			}
		}
		
		// this loop checks that the leading coefficients are 1 AND that they have nothing above them,
		//  all in one fell swoop, by checking that all the pivot columns are unit vectors.
		for(int row=1; row<height; row++) {
			if( pivotCols[row] == width) {
				// success because we've gotten to the end and had no failures.
				// I have to end it early here to prevent an IndexOutOfBoundsException
				return true;
			}
			if(! getCol(pivotCols[row]).isUnitVector()  ) {
				// failed because either there is not a 1, or there is stuff other than the 1.
				return false;
			}
		}
		// If there are no zero rows the code could reach here. No failures means success!
		return true;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// Time for the actual refing and rrefing.
	
	public static Matrix ref(Matrix M) {
		int height = M.height;
		int width = M.width;
		
		// we will keep looking at smaller and smaller parts of the matrix.
		int firstRow = 0;
		int firstCol = 0;
		
		while(firstRow < height && firstCol < width) {
		
			// we need to find the leftmost nonzero column and then ignore all columns to the left of that
			int pivotCol = width;
			for(int col=firstCol; col<width; col++) {
				if( M.getCol(col).nonzeroPast(firstRow) ) {
					pivotCol = col;
					break;
				}
			}
			firstCol = pivotCol;
			
			if(firstCol == width) {
				// the matrix is entirely zeroes, and is already in ref.
				break;
			}
			
			// now we need to find the first nonzero entry of this column, put it at the top, and make it 1.
			int firstNonzeroEntry = height;
			// ^this will get overwritten, but I have to give it an initial value
			//  so that the compiler doesn't yell at me.
			for (int row=firstRow; row<height; row++) {
				if(M.get(row,firstCol).nonzero()) {
					firstNonzeroEntry = row;
					break;
				}
			}
			M = M.rowSwap(firstRow, firstNonzeroEntry); // if these are both 0, nothing happens
			M = M.rowScale(firstRow, M.get(firstRow,firstCol).inverse());
			M = M.cleanup();
			
			// now we need to make sure there are zeroes below the pivot entry.
			for(int row = firstRow+1; row<height; row++) {
				M = M.rowAdd(row, firstRow, M.get(row,firstCol).negate());
				M = M.cleanup();
			}
			
			// now the top row is good. move down one row and repeat.
			firstRow++;
		}
		// the matrix will now be ref with the pivot entries all 1.
		return M.cleanup();
	}
	
	// a version that prints out the matrix after each step
	public static Matrix refDebug(Matrix M) {
		int height = M.height; int width = M.width;
		int firstRow = 0; int firstCol = 0;
		while(firstRow < height && firstCol < width) {
			int pivotCol = width;
			for(int col=firstCol; col<width; col++)
				if( M.getCol(col).nonzeroPast(firstRow) ) {pivotCol = col; break;}
			firstCol = pivotCol;
			if(firstCol == width) break;
			int firstNonzeroEntry = height;
			for (int row=firstRow; row<height; row++)
				if(M.get(row,firstCol).nonzero()) {firstNonzeroEntry = row; break;}
			
			M = M.rowSwap(firstRow, firstNonzeroEntry);
			M = M.rowScale(firstRow, M.get(firstRow,firstCol).inverse());
			M = M.cleanup();
			System.out.println(M);
			
			for(int row = firstRow+1; row<height; row++) {
				M = M.rowAdd(row, firstRow, M.get(row,firstCol).negate());
				M = M.cleanup();
				System.out.println(M);
			}
			
			firstRow++;
		}
		return M.cleanup();
	}
	
	
	public static Matrix rref(Matrix M) {
		M = ref(M);
		
		int[] pivotCols = new int[M.height];
		for(int i=0;i<M.height;i++) {
			pivotCols[i] = M.pivotCol(i);
		}
		
		// go backwards through the pivot entries and make sure each one has zeroes above it
		for(int i=M.height-1; i>=0; i--) {
			if(pivotCols[i] == M.width) {
				// this indicates that the row is empty.
				continue;
			}
			
			// subtract from each row above this one.
			for(int j=i-1; j>=0; j--) {
				M = M.rowAdd(j, i, M.get(j,pivotCols[i]).negate());
			}
			
		}
		// that's all!
		// thanks for coming!
		return M.cleanup();
	}
	
	// non-static versions
	public Matrix ref() {
		return Matrix.ref(this);
	}
	public Matrix rref() {
		return Matrix.rref(this);
	}
	
	
	
	
	
	
	
	
	// One thing you can do with rref is solve Mv = u for vectors M and u. A special case of this is
	//  finding the kernel (null space) of a matrix.
	
	// dim of the row space
	public int rank() {
		int rank = 0;
		for(Vector row: rref(this).rows()) {
			if(row.nonzero()) {
				rank++;
			}
		}
		return rank;
	}
	
	// dim of the kernel
	public int nullity() {
		return width - rank();
		// by the rank-nullity theorem
	}
	
	// static versions
	public static int rank(Matrix M) {
		return M.rank();
	}
	public static int nullity(Matrix M) {
		return M.nullity();
	}
	
	// finds a basis for the kernel.
	// this method is static so that I don't accidentally reference entries of the matrix when I mean
	//  to mention entries of rref(M).
	public static Vector[] kernel(Matrix M) {
		// each row in the matrix can be interpreted as an equation:
		//  [0 0 ... 1 X Y Z] means
		//  1*w + X*x + Y*y + Z*z = 0
		//  so if we know x, y, and z, we can solve for w.
		
		M = rref(M);
		
		int nullity = M.nullity();
		if(nullity==0) {
			return new Vector[0];
		}
		
		// We need to find which columns are never pivot columns.
		// These columns correspond to free variables.
		int[] freeVars = new int[nullity];
		for(int col=0, k=0, row=0; col<M.width; col++) {
			if(row >= M.height) {
				freeVars[k] = col;
				k++;
			} else if(M.pivotCol(row) == col) {
				row++;
			} else {
				freeVars[k] = col;
				k++;
			}
		}
		
		
		Vector[] kernel = new Vector[nullity];
		// Now it's time to find the actual vectors.
		for(int k=0; k<nullity; k++) {
			Complex[] vec = new Complex[M.width];
			
			for(int l=0; l<nullity; l++) {
				// for each vector, set one of the free variables to one and the others to zero.
				vec[freeVars[l]] = (l==k) ? Complex.ONE: Complex.ZERO;
			}
			
			for(int row=M.height-1; row>=0; row--) {
				
				if(!M.getRow(row).nonzero()) {
					continue;
					// rows that are all zero give us no information.
				}

				int index = M.pivotCol(row);
				vec[index] = Complex.ZERO;
				for(int i=index+1; i<M.width; i++) {
					vec[index] = Complex.add(vec[index], Complex.mul(M.get(row,i),vec[i]) );
					// 1*w + X*x + Y*y + Z*z = 0
					// w = - (X*x + Y*y + Z*z)
					// w = - Σ(get(row,i)*vec[i])
				}
				vec[index] = vec[index].negate();
			}
			
			kernel[k] = new Vector(vec);
		}
		
		return kernel;
	}
	
	// non-static version
	public Vector[] kernel() {
		return kernel(this);
	}
	
	// returns one vector from the kernel.
	public Vector kernelExample() {
		Vector[] ker = kernel();
		if(ker.length == 0) {
			return Vector.ZERO(width);
		} else {
			return ker[0];
		}
	}
	
	// returns one vector from the preimage of target.
	// this function works very similarly to kernel.
	// if target is not in the range of the matrix, this method returns null.
	public Vector preimage(Vector target) {
		
		if(target.size != this.height) {
			throw new DimensionMismatchException(
					"The given vector is not in the codomain of the matrix. The heights must match.");
		}
		
		// Av₂=v₁   =>    v₂=A⁻¹v₁
		if(isInvertible()) {
			return inverse().mul(target);
		}
		
		// otherwise, construct the augmented matrix [ A | v₁ ]
		Vector[] augmentedCols = new Vector[width+1];
		for(int i=0; i<width; i++) {
			augmentedCols[i] = getCol(i);
		}
		augmentedCols[width] = target;
		Matrix augmentedMatrix = Matrix.fromCols(augmentedCols);
		augmentedMatrix = rref(augmentedMatrix);
		
		int nullity = augmentedMatrix.nullity();
		
		// We need to find which columns are never pivot columns.
		// These columns correspond to free variables.
		int[] freeVars = new int[nullity];
		for(int col=0, k=0, row=0; col<width; col++) {
			if(row >= height) {
				freeVars[k] = col;
				k++;
			} else if(augmentedMatrix.pivotCol(row) == col) {
				row++;
			} else {
				freeVars[k] = col;
				k++;
			}
		}
		
		Complex[] vec = new Complex[width];
		
		for(int l=0; l<nullity; l++) {
			// set all of the free variables to zero to make things easier.
			// the vector will still end up nonzero unless target is zero.
			vec[freeVars[l]] = Complex.ZERO;
		}
		
		for(int row=height-1; row>=0; row--) {
			
			if(!augmentedMatrix.getRow(row).nonzero()) {
				continue;
				// rows that are all zero give us no information.
			}

			int index = augmentedMatrix.pivotCol(row);
			if(index == width) {
				// this means that one of the rows is [0, 0, ... 0 | V] for nonzero V
				// 0 + 0 + ... + 0 = V is a contradiction
				// returning null is better than raising an error because this information can be useful
				//  and we shouldn't expect the user to know whether a particular vector is in the range
				//  of a matrix.
				return null;
			}
			vec[index] = augmentedMatrix.get(row, width);
			for(int i=index+1; i<width; i++) {
				vec[index] = Complex.sub(vec[index], Complex.mul(augmentedMatrix.get(row,i),vec[i]) );
				// 1*w + X*x + Y*y + Z*z = V
				// w = V - (X*x + Y*y + Z*z)
				// w = get(row,width) - Σ(get(row,i)*vec[i])
			}
		}
		
		return new Vector(vec);
		
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// necessary to find the determinant using the minors method.
	public Matrix minor(int strikeRow, int strikeCol) {
		Complex[][] newArray = new Complex[height-1][width-1];
		
		for(int row=0; row<strikeRow; row++) {
			for(int col=0; col<strikeCol; col++) {
				newArray[row][col] = array[row][col];
			}
			for(int col=strikeCol+1; col<width; col++) {
				newArray[row][col-1] = array[row][col];
			}
		}
		for(int row=strikeRow+1; row<height; row++) {
			for(int col=0; col<strikeCol; col++) {
				newArray[row-1][col] = array[row][col];
			}
			for(int col=strikeCol+1; col<width; col++) {
				newArray[row-1][col-1] = array[row][col];
			}
		}
		return new Matrix(newArray);
	}
	
	
	
	// This method finds the determinant using expansion by minors.
	public Complex determinant() {
		if(!isSquare()) {
			throw new DimensionMismatchException("Only square matrices have determinants.");
		}
		
		// base case since this definition is recursive
		if(width == 1) {
			return get(0,0);
		}
		
		Complex sum = Complex.ZERO;
		int sign = 1;
		for(int i=0; i<width; i++) {
			// sum +=   M[0][i] * det(M_0,i) * (-1)^i
			sum = Complex.add(sum, Complex.mul(get(0,i), det(minor(0,i))).scale(sign) );
			
			sign *= -1;
		} 
		
		return sum;
	}
	// static version
	public static Complex det(Matrix M) {
		return M.determinant();
	}
	
	public boolean isInvertible() {
		if(!isSquare()) {
			return true;
		}
		return determinant().nonzero();
	}
	
	
	// Now that rref is finished, we can find the inverse using the row reduction method.
	// Keep track of two matrices. One starts out as M and the other starts out as the identity.
	//  Row reduce M until it becomes the identity, and each operation you do to M, also do to the
	//  second matrix. By the end, the second matrix will be the identity.
	// This is a static method so that I don't accidentally reference "this". The non-static version
	//  will just call the static one.
	public static Matrix invert(Matrix M) {
		if(!M.isSquare()) {
			throw new DimensionMismatchException("Only square matrices can have inverses.");
		}
		if(!M.isInvertible()) {
			throw new ArithmeticException("Attempt to invert a non-invertible matrix.");
		}
		
		int height = M.height;
		int width = M.width;
		
		// here is the second matrix, that will eventually become the inverse.
		Matrix I = IDENTITY(M.width);
		
		// now for the repeat of ref() and rref()
		
		int firstRow = 0;
		int firstCol = 0;
		
		while(firstRow < height && firstCol < width) {
			int pivotCol = width;
			for(int col=firstCol; col<width; col++) {
				if( M.getCol(col).nonzeroPast(firstRow) ) {
					pivotCol = col;
					break;
				}
			}
			firstCol = pivotCol;
			
			if(firstCol == width) break;
			
			int firstNonzeroEntry = height;
			for (int row=firstRow; row<height; row++) {
				if(M.get(row,firstCol).nonzero()) {
					firstNonzeroEntry = row;
					break;
				}
			}
			// here, I do the same thing to the second matrix.
			I = I.rowSwap(firstRow, firstNonzeroEntry);
			M = M.rowSwap(firstRow, firstNonzeroEntry);
			I = I.rowScale(firstRow, M.get(firstRow,firstCol).inverse());
			M = M.rowScale(firstRow, M.get(firstRow,firstCol).inverse());
			
			for(int row = firstRow+1; row<height; row++) {
				I = I.rowAdd(row, firstRow, M.get(row,firstCol).negate());
				M = M.rowAdd(row, firstRow, M.get(row,firstCol).negate());
			}
			firstRow++;
		}
		M = M.cleanup();
		I = I.cleanup();
		
		// that was ref, now here's the rest of rref.
		
		int[] pivotCols = new int[M.height];
		for(int i=0;i<M.height;i++) {
			pivotCols[i] = M.pivotCol(i);
		}
		
		for(int i=M.height-1; i>=0; i--) {
			if(pivotCols[i] == M.width)	continue;
			
			for(int j=i-1; j>=0; j--) {
				I = I.rowAdd(j, i, M.get(j,pivotCols[i]).negate());
				M = M.rowAdd(j, i, M.get(j,pivotCols[i]).negate());
			}
			
		}
		
		// and we're done!
		return I.cleanup();
		
	}
	
	// here's the non-static version.
	public Matrix inverse() {
		return invert(this);
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// eigenpalooza
	
	public Polynomial charpoly() {
		if(!isSquare()) {
			throw new DimensionMismatchException("Only square matrices have characteristic polynomials.");
		}
		
		PolynomialMatrix AMinusLambdaI = PolynomialMatrix.sub( new PolynomialMatrix(this),
				PolynomialMatrix.IDENTITY(this.height).scale(Polynomial.X) );
		return AMinusLambdaI.determinant();
	}
	public static Polynomial charpoly(Matrix M) {
		return M.charpoly();
	}
	
	public Complex[] eigenvalues() {
		if(!isSquare()) {
			throw new DimensionMismatchException("Only square matrices have eigenpairs.");
		}
		
		return Algebra.solve(charpoly());
		
	}
	
	// returns one of the eigenvectors associated with the eigenvalue lambda.
	public Vector eigenvector(Complex lambda) {
		if(!isSquare()) {
			throw new DimensionMismatchException("Only square matrices have eigenpairs.");
		}
		// we need to find the vector v such that
		// Mv = λv
		// (M-λI)v = 0
		// so we need to find a vector in the null space of (M-λI)
		
		Matrix R = rref( sub(this, IDENTITY(height).scale(lambda)) );
		if(R.equals(IDENTITY(height))) {
			// this is not an eigenvector
			throw new IllegalArgumentException(lambda + " is not an eigenvector of the matrix.");
		}
		
		Complex[] evec = new Complex[this.height];
		evec[height-1] = Complex.ONE;
		for(int i=this.height-2; i>=0; i--) {
			// if the ith row is [0,0...1,W,X,Y,Z] then the ith entry of the eigenvector should be 
			// -(W*w + X*x + Y*y + Z*z)
			evec[i] = Complex.ZERO;
			for(int j=i+1; j<height; j++) {
				evec[i] = Complex.add( evec[i] , Complex.mul(evec[j], R.get(i, j))  );
			}
			evec[i] = evec[i].negate();
			// now evec[i] is the correct sum
		}
		
		return new Vector(evec);
	}
	
	// returns a basis for the subspace of eigenvectors associated with the eigenvalue lambda
	public Vector[] eigenvectors(Complex lambda) {
		if(!isSquare()) {
			throw new DimensionMismatchException("Only square matrices have eigenpairs.");
		}
		
		// ker(A-λI)
		return kernel( sub(this, IDENTITY(height).scale(lambda))  );
	}
	
	
	
	
	
	
	
	
	
	
	public Matrix[] diagonalize() {
		if(!isSquare()) {
			throw new DimensionMismatchException("Only square matrices can be diagonalized.");
		}
		
		Complex[] lambda = eigenvalues();
		
		// get rid of repeated eigenvalues
		for(int i=0; i<height-1; i++) {
			for(int j=i+1; j<height; j++) {
				if(lambda[i]==null || lambda[j]==null) {
					continue;
				}
				if(lambda[j].equals(lambda[i])) {
					lambda[j] = null;
				}
			}
		}

		// now we need to make two lists, of e.vecs and their respective e.vals .
		Complex[] evalsInOrder = new Complex[height];
		Vector[] evecs = new Vector[height];
		int numEvecs = 0;
		
		for(int i=0; i<height; i++) {
			if(lambda[i] == null) {
				continue;
			}
			
			for(Vector v: eigenvectors(lambda[i])) {
				evecs[numEvecs] = v;
				evalsInOrder[numEvecs] = lambda[i];
				numEvecs++;
			}
		}
		
		// if the matrix is not diagonalizable, we will end up with too few e.vecs .
		if(numEvecs != height) {
			throw new ArithmeticException("The matrix is not diagonalizable.");
		}
		
		Matrix P = Matrix.fromCols(evecs);
		Matrix D = diag(evalsInOrder);
		
		return new Matrix[] {P, D, P.inverse()};
	}
	
	
	
	public Matrix[] jordanForm() {
		if(!isSquare()) {
			throw new DimensionMismatchException("Only square matrices have Jordan forms.");
		}
		
		Complex[] lambda = eigenvalues();
		int[] multiplicity = new int[height]; 
		// the array is filled with 0's by default. it should be ones.
		for(int i=0; i<height; i++) multiplicity[i] = 1;
		
		// get rid of repeated eigenvalues, and keep track of the multiplicity of each
		for(int i=0; i<height-1; i++) {
			for(int j=i+1; j<height; j++) {
				if(lambda[i]==null || lambda[j]==null) {
					continue;
				}
				if(lambda[j].equals(lambda[i])) {
					lambda[j] = null;
					multiplicity[i]++;
				}
			}
		}
		
		// now we will start putting together J (the jordan-block matrix) and the columns of P.
		Complex[][] JArray = new Complex[height][height];
		Vector[] vecs = new Vector[height];
		int numVecs = 0;
		
		// loop through eigenvalues
		for(int i=0; i<height; i++) {
			if(lambda[i] == null) {
				continue;
			}
			
			Matrix AMinusLambdaI = sub(this, IDENTITY(height).scale(lambda[i]));
			
			for(Vector v: eigenvectors(lambda[i])) {
				JArray[numVecs][numVecs] = lambda[i];
				vecs[numVecs] = v;
				numVecs++;
				
				// try generalized eigenvectors
				v = AMinusLambdaI.preimage(v);
				while(v != null) {
					JArray[numVecs][numVecs] = lambda[i];
					JArray[numVecs-1][numVecs] = Complex.ONE;
					vecs[numVecs] = v;
					numVecs++;
					v = AMinusLambdaI.preimage(v);
				}
				
			}
		}
		
		// fill in the null values
		for(int row=0;row<height;row++) {
			for(int col=0;col<height;col++) {
				if(JArray[row][col] == null) {
					JArray[row][col] = Complex.ZERO;
				}
			}
		}
		
		Matrix P = Matrix.fromCols(vecs);
		Matrix J = new Matrix(JArray);
		
		return new Matrix[] {P, J, P.inverse()};
	}
	
	
	public Matrix jordanCanonicalForm() {
		
		// this algorithm relies on an implementation detail: the J returned by jordanForm() always
		//  puts blocks with the same eigenvalue next to one another.
		
		if (height == 0 && width == 0) {
			return Matrix.ZERO(0);
		}
		
		Matrix J = jordanForm()[1]; // square-ness is checked here
		Map<Complex, ArrayList<Integer>> blocks = new TreeMap<Complex, ArrayList<Integer>>();
		
		Complex eval = J.get(0, 0);
		blocks.put(eval, new ArrayList<>());
		blocks.get(eval).add(1);
		
		for(int i=1; i<height; i++) {
			eval = J.get(i, i);
			
			if(! eval.equals(J.get(i-1, i-1))) {
				blocks.put(eval, new ArrayList<>());
				blocks.get(eval).add(1);
				continue;
			}
			// if we get here, this e.val is the same as the last one.
			ArrayList<Integer> blockSizes = blocks.get(eval);
			if(J.get(i-1, i).equals(Complex.ONE)) {
				int numBlocks = blockSizes.size();
				blockSizes.set(numBlocks, blockSizes.get(numBlocks) + 1);
			} else {
				// this is the start of a new block with the same eigenvalue
				blockSizes.add(1);
			}
		}
		
		// now we have this treemap with eigenvalue: {sizes of jordan blocks}.
		// to ensure that this is invariant, we need to sort the lists.
		for(Complex c: blocks.keySet()) {
			blocks.get(c).sort(null);
		}
		
		// now construct the actual jordan blocks
		int numJordanBlocks = J.height - (int) J.minor(height-1, 0).trace().re;
		Matrix[] blocksAsMatrices = new Matrix[numJordanBlocks];
		int progress = 0;
		for(Complex c: blocks.keySet()) {
			for(int i: blocks.get(c)) {
				blocksAsMatrices[progress] = jordanBlock(i, c);
				progress++;
			}
		}
		
		return diag(blocksAsMatrices);	
	}
	
	
	public boolean isSimilar(Matrix that) {
		return this.jordanCanonicalForm().equals(that.jordanCanonicalForm());
	}
	
	
	
	
	
	
	
	// the diagonal and nilpotent part of a matrix in jordan-block form
	private Matrix diagonalPartJ() {
		int size = height;
		Complex[][] array = new Complex[size][size];
		for(int row=0; row<size; row++) {
			for(int col=0; col<size; col++) {
				array[row][col] = (row==col) ? get(row,col) : Complex.ZERO;
			}
		}
		return new Matrix(array);
	}
	private Matrix nilpotentPartJ() {
		int size = height;
		Complex[][] array = new Complex[size][size];
		for(int row=0; row<size; row++) {
			for(int col=0; col<size; col++) {
				array[row][col] = (col==row+1) ? get(row,col) : Complex.ZERO;
			}
		}
		return new Matrix(array);
	}
	
	// the diagonal and nilpotent part of a generic square matrix.
	// the matrices D and N such that A=D+N and DN=ND
	public Matrix diagonalPart() {
		if(!isSquare()) {
			throw new DimensionMismatchException(
					"Only square matrices can be broken into diagonal and nilpotent parts.");
		}
		Matrix[] jf = jordanForm();
		return mul(jf[0], jf[1].diagonalPartJ(), jf[2]);
	}
	public Matrix nilpotentPart() {
		if(!isSquare()) {
			throw new DimensionMismatchException(
					"Only square matrices can be broken into diagonal and nilpotent parts.");
		}
		Matrix[] jf = jordanForm();
		return mul(jf[0], jf[1].nilpotentPartJ(), jf[2]);
	}
	
	
	
	
	
	// exp of a jordan-block matrix
	private static Matrix expJ(Matrix J) {
		if(!J.isSquare()) {
			throw new DimensionMismatchException(
					"Only square matrices can be exponentiated.");
		}
		
		int size = J.height;
		
		Matrix D = J.diagonalPartJ();
		Matrix N = J.nilpotentPartJ();
		// exp(J) = exp(D) * Σ_(k=0)^(n-1) (N^k)/k!
		
		Matrix expN = ZERO(size);
		for(int k=0; k<size; k++) {
			expN = add(expN, N.pow(k).scale(1.0/factorial(k)) );
		}
		
		Complex[][] expD = new Complex[size][size];
		for(int row=0; row<size; row++) {
			for(int col=0; col<size; col++) {
				expD[row][col] = (row==col) ? Complex.exp(D.get(row,col)) : Complex.ZERO;
			}
		}
		
		return mul(new Matrix(expD),expN);
	}
	
	// exp of any square matrix
	public static Matrix exp(Matrix M) {
		if(!M.isSquare()) {
			throw new DimensionMismatchException(
					"Only square matrices can be exponentiated.");
		}
		Matrix[] jf = M.jordanForm();
		return mul(jf[0], expJ(jf[1]), jf[2]);
	}
	
	
	// sin and cos of a jordan-block matrix
	private static Matrix sinJ(Matrix J) {
		if(!J.isSquare()) {
			throw new DimensionMismatchException(
					"Only square matrices can be plugged into functions.");
		}
		
		int size = J.height;
		
		Matrix D = J.diagonalPartJ();
		Matrix N = J.nilpotentPartJ();
		// exp(J) = Σ_(k=0)^(n-1) [sin⁽ᵏ⁾(D)(N^k)/k!]
		
		Matrix sum = ZERO(size);
		for(int k=0; k<size; k++) {
			sum = add(sum, mul(sinDerivDiag(k,D), N.pow(k).scale(1.0/factorial(k))) );
		}
		return sum;
	}
	private static Matrix cosJ(Matrix J) {
		if(!J.isSquare()) {
			throw new DimensionMismatchException(
					"Only square matrices can be plugged into functions.");
		}
		
		int size = J.height;
		
		Matrix D = J.diagonalPartJ();
		Matrix N = J.nilpotentPartJ();
		// exp(J) = Σ_(k=0)^(n-1) [sin⁽ᵏ⁾(D)(N^k)/k!]
		
		Matrix sum = ZERO(size);
		for(int k=0; k<size; k++) {
			sum = add(sum, mul(sinDerivDiag(k+1,D), N.pow(k).scale(1.0/factorial(k))) );
		}
		return sum;
	}
	
	// sin and cos of any square matrix
	public static Matrix sin(Matrix M) {
		if(!M.isSquare()) {
			throw new DimensionMismatchException(
					"Only square matrices can be plugged into functions.");
		}
		Matrix[] jf = M.jordanForm();
		return mul(jf[0], sinJ(jf[1]), jf[2]);
	}
	public static Matrix cos(Matrix M) {
		if(!M.isSquare()) {
			throw new DimensionMismatchException(
					"Only square matrices can be plugged into functions.");
		}
		Matrix[] jf = M.jordanForm();
		return mul(jf[0], cosJ(jf[1]), jf[2]);
	}
	
	
	private static Matrix arcsinJ(Matrix J) {
		if(!J.isSquare()) {
			throw new DimensionMismatchException(
					"Only square matrices can be plugged into functions.");
		}
		
		int size = J.height;
		
		Matrix D = J.diagonalPartJ();
		Matrix N = J.nilpotentPartJ();
		// exp(J) = Σ_(k=0)^(n-1) [sin⁽ᵏ⁾(D)(N^k)/k!]
		
		Matrix sum = ZERO(size);
		for(int k=0; k<size; k++) {
			sum = add(sum, mul(arcsinDerivDiag(k,D), N.pow(k).scale(1.0/factorial(k))) );
		}
		return sum;
	}
	public static Matrix arcsin(Matrix M) {
		if(!M.isSquare()) {
			throw new DimensionMismatchException(
					"Only square matrices can be plugged into functions.");
		}
		Matrix[] jf = M.jordanForm();
		return mul(jf[0], arcsinJ(jf[1]), jf[2]);
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// useful
	// from https://stackoverflow.com/questions/388461/how-can-i-pad-a-string-in-java
	private static String padRight(String s, int n) {
	     return String.format("%-" + n + "s", s);  
	}
	
	
	// sin⁽ᵏ⁾ of a diagonal matrix D
	private static Matrix sinDerivDiag(int k, Matrix D) {
		Complex[][] array = new Complex[D.height][D.height];
		for(int row=0; row<D.height; row++) {
			for(int col=0; col<D.height; col++) {
				array[row][col] = (row==col) ? Complex.sinDeriv(k, D.get(row,col)) : Complex.ZERO;
			}
		}
		return new Matrix(array);
	}
	
	private static Matrix arcsinDerivDiag(int k, Matrix D) {
		Complex[][] array = new Complex[D.height][D.height];
		for(int row=0; row<D.height; row++) {
			for(int col=0; col<D.height; col++) {
				array[row][col] = (row==col) ? Analysis.arcsinDeriv(k, D.get(row,col)) : Complex.ZERO;
			}
		}
		return new Matrix(array);
	}
	
	
	
	
}
