package ccd.model.BCCD;

public class Utils {
    public static double logOrZero(double value) {
        return value == 0.0 ? 0.0 : Math.log(value);
    }
}
