package linalg;


public class Main {
	@SuppressWarnings("unused")
	public static void main(String[] args) {
		
		int[][] array = new int[][] {
			{2,3,4},
			{1,-2,0},
			{132,-5,10},
			{2,-4,0}
		};
		
		int[][] array1 = new int[][] {
			{-5,4,1,-1},
			{-3,2,0,2},
			{4,0,4,-4},
			{-2,2,1,0}
		};
		
		int[][] array2 = new int[][] {
			{2,1,0,-3},
			{1,1,1,-3},
			{1,0,1,-2},
			{1,0,-1,0}
		};
		
		int[][] array3 = new int[][] {
			{2,4,-1,2},
			{4,0,10,5},
			{-4,1,1,-1},
			{1,5,-2,1}
		};
		
		
		double[][] array4 = new double[][] {
			{0,1,1,1,0,1,0},
			{1,0,1,0,1,0,0},
			{1,1,0,0,0,0,1},
			{1,0,0,0,0,0,0},
			{0,1,0,0,0,0,0},
			{1,0,0,0,0,0,0},
			{0,0,1,0,0,0,0}
		};
		
		
		Matrix M;
		Matrix M1 = new Matrix(array1);
		Matrix M2 = new Matrix(array2);
		Matrix M3 = new Matrix(array3);
		Matrix M4 = new Matrix(array4);
		
		
		
		System.out.println(M4.jordanForm());
		
	}
}
