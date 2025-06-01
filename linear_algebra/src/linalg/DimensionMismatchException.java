package linalg;

public class DimensionMismatchException extends IllegalArgumentException {

	public DimensionMismatchException(String string) {
		super(string);
	}
	
	public DimensionMismatchException() {
		super("Dimensions of matrix different than expected.");
	}

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

}
