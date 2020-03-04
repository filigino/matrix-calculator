/*
 * Name: Gino Salayo
 * Description: This program returns the inverse of the product of a set of matrices.
 * Input Format (space-separated integers):
 * N (the number of matrices)
 * m1, n1 (size of first matrix)
 * m2, n2 (size of second matrix)
 * ...
 * mN, nN (size of last matrix)
 * m1 x n1 elements for first matrix row by row
 * m2 x n2 elements for second matrix row by row
 * ...
 * mN x nN elements for nth matrix row by row
 */

public class MatrixCalculator {
	// Computes cofactor matrix with determinant of submatrix obtained by
	// omitting ith row and jth column
	public static double[][] findCofactors(double[][] matrix) {
		double[][] cofactors = new double[matrix.length][matrix.length];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix.length; j++) {
				cofactors[i][j] = Math.pow(-1, i + j)
				    *findDet(findSubmatrix(matrix, i, j));
			}
		}
		return cofactors;
	}
	
	public static double findDet(double[][] matrix) {
		int l = matrix.length;
		if (l == 2) {
			double det = (matrix[0][0])*(matrix[1][1]) - (matrix[0][1])*(matrix[1][0]);
				//calculate determinant of 2x2 matrix with a*d - b*c and store in double det
			return det;
		}
		else { //if matrix is larger than 2x2
			double det = 0; //initialize det = 0 to store accumulated value of determinant
			for (int i = 0; i < matrix.length; i++) { //for loop iterates for each column of input matrix
				det += Math.pow(-1,i)*matrix[0][i]*findDet(findSubmatrix(matrix,0,i));
					//adds to running value of determinant by multiplying (-1)^n with the ith entry of the first
					//row of the input matrix and with the determinant of the submatrix obtained by omitting the
					//1st row and ith column
			}
			return det;
		}
	}

	public static double[][] findInverse(double[][] matrix) { //calculates determinant of input matrix
		double[][] inverse = new double[matrix.length][matrix.length];
			//creates empty matrix with same dimensions as input matrix
		double det = findDet(matrix); //finds determinant of input matrix
		
		if (matrix.length == 1) {
			inverse[0][0] = 1/matrix[0][0]; //calculates inverse of 1x1 matrix by taking inverse of the entry
		}
		else if (matrix.length == 2) {
			inverse[0][0] = matrix[1][1]/det + 0.0;
			inverse[0][1] = -matrix[0][1]/det + 0.0;
			inverse[1][0] = -matrix[1][0]/det + 0.0;
			inverse[1][1] = matrix[0][0]/det + 0.0;
				//calculates inverse of 2x2 matrix with 1/det*([d, -b], [-c, a])
		}
		else {
			double[][] cofactor = findCofactors(matrix); //finds cofactor of input matrix
			for (int i = 0; i < inverse.length; i++) { //for loop iterates for each row of input matrix
				for (int j = 0; j < inverse.length; j++) { //for loop iterates for each column of input matrix
					inverse[i][j] = cofactor[j][i]/det + 0.0;
						//calculates i, jth entry by calculating adjoint of input matrix and dividing by determinant
				}
			}
		}
		return inverse;
	}
	
	public static double[][] findSubmatrix(double[][] matrix, int m, int n) {
		//returns submatrix obtained by omitting mth row and nth column of input matrix
		double[][] submatrix = new double[matrix.length - 1][matrix.length - 1];
			//creates empty matrix of size (m-1)x(n-1)
		
		for (int j = 0; j < m; j++) { //for loop iterates for each row before omitted row
			for (int k = 0; k < n; k++) { //for loop iterates for each column before omitted column
				submatrix[j][k] = matrix[j][k]; //stores submatrix values
			} //end for loop
			for (int k = submatrix.length; k > n; k--) {
				//for loop iterates for each column after omitted column, starting from right side of input matrix
				submatrix[j][k-1] = matrix[j][k]; //stores submatrix values
			} //end for loop
		} //end for loop
		
		for (int j = submatrix.length; j > m; j--) {
			//for loop iterates for each row after omitted row, starting from the bottom of input matrix
			for (int k = 0; k < n; k++) { //for loop iterates for each column before omitted column
				submatrix[j-1][k] = matrix[j][k]; //stores submatrix values
			} //end for loop
			for (int k = submatrix.length; k > n; k--) {
				//for loop iterates for each column after omitted column, starting from right side of input matrix
				submatrix[j-1][k-1] = matrix[j][k]; //stores submatrix values
			} //end for loop
		} //end for loop
		return submatrix; //returns submatrix
	} //end findSubmatrix
	
	public static double[][] multiplyMatrices(double[][] a, double[][] b) { 
		//multiplies two matrices together
		int m1 = a.length; //stores number of rows of first matrix in m1
		int n1 = a[0].length; //stores number of columns of first matrix in n1
		int n2 = b[0].length; //stores number of rows of second matrix in n2
		double[][] c = new double[m1][n2]; //creates empty matrix of size m1xn2

		for (int i = 0; i < m1; i++) { //for loop iterates for each row of first matrix
			for (int j = 0; j < n2; j++) { //for loop iterates for each column of second matrix
				for (int k = 0; k < n1; k++) { //for loop iterates for each column of first matrix
					c[i][j] += a[i][k] * b[k][j];
						//finds i, jth entry of product of the matrices by taking dot product of ith row of first
						//matrix and jth column of second matrix
				} //end for loop
			} //end for loop
		} //end for loop
		return c; //returns product of matrices
	} //end multiplyMatrices
	
	public static void main(String[] args) {
		int N = Integer.parseInt(args[0]); // number of matrices
		
		for (int i = 0; i < (N - 1); i++) {
			//sanity check to ensure adjacent matrices have matching dimensions to be multiplied together
			int a = Integer.parseInt(args[2 + 2*i]); //stores number of columns of matrix in a
			int b = Integer.parseInt(args[3 + 2*i]); //stores number of rows of next matrix in b
			
			if (a != b) {
				//if number of columns of a matrix does not match number of rows of adjacent matrix,
				//matrices cannot be multiplied together
				System.out.println("Multiplication error."); //prints error message
				return; //if two matrices cannot be multiplied together, exits the entire program
			} //end if statement
		} //end for loop
		
		double[][][] matrices = new double[N][][]; //creates 3d array of doubles to store all input matrices
		int counter = 2*N + 1; //creates counter to keep track of location in array args
		
		for (int i = 0; i < N; i++) { //for loop iterates for each input matrix
			int m = Integer.parseInt(args[2*i + 1]);
			int n = Integer.parseInt(args[2*i + 2]);
				//stores number of rows and columns of ith matrix in m and n
			matrices[i] = new double[m][n]; //creates 2d array for matrix in ith entry of matrices
			
			for (int j = 0; j < m; j++) { //for loop iterates for each row of ith matrix
				for (int k = 0; k < n; k++) { //for loop iterates for each column of ith matrix
					matrices[i][j][k] = Double.parseDouble(args[counter]);
					counter++; //increments counter by 1
						//stores entries of the matrices into corresponding arrays for each input matrix
				} //end for loop
			} //end for loop
		} //end for loop
		
		double[][] product = matrices[0]; //creates new matrix for product of two matrices
		
		for (int i = 1; i < N; i++) { //for loop iterates for each input matrix
			double[][] a = matrices[i] ; //stores ith matrix in a
			double[][] b = multiplyMatrices(product, a); //multiplies adjacent matrices
			product = b; //stores product in 2d array product
		}
		
		if (product.length != product[0].length) { //if product is not a square matrix
			System.out.println("Matrix not invertible.");
			return; //stops program
		}
		else if (findDet(product) == 0 && product.length != 1) { //if determinant is 0 and product is not 1x1
			System.out.println("Matrix not invertible.");
			return; //stops program
		}
		
		double[][] inverse = findInverse(product);

		for (int i = 0; i < inverse.length; i++) { //for loop iterates for each row of inverse
			for (int j = 0; j < inverse.length; j++) { //for loop iterates for each colummn of inverse
				System.out.print(inverse[i][j] + " "); //prints each element of inverse
			}
		}		
	}
}
