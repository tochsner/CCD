package ccd.tools;

import org.apache.commons.math3.linear.RealMatrix;

public class MatrixUtils {
    /**
     * Fills the given array with the values in the given matrix.
     */
    public static void fillArray(RealMatrix matrix, double[][] array) {
        for (int i = 0; i < matrix.getRowDimension(); i++) {
            for (int j = 0; j < matrix.getColumnDimension(); j++) {
                array[i][j] = matrix.getEntry(i, j);
            }
        }
    }
}
