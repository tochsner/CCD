package ccd.model;

import beast.base.evolution.tree.Node;
import javafx.scene.input.ScrollEvent;

public class Utils {
    public static double logOrZero(double value) {
        return value == 0.0 ? 0.0 : Math.log(value);
    }

    public static <T extends CladePartition> Clade getFirstClade(T partition) {
        return Utils.getFirstClade(partition.getChildClades()[0], partition.getChildClades()[1]);
    }

    public static <T extends Clade> T getFirstClade(T clade1, T clade2) {
        if (clade1.getCladeInBits().nextSetBit(0) < clade2.getCladeInBits().nextSetBit(0)) {
            return clade1;
        } else {
            return clade2;
        }
    }

    public static <T extends CladePartition> Clade getSecondClade(T partition) {
        return Utils.getSecondClade(partition.getChildClades()[0], partition.getChildClades()[1]);
    }

    public static <T extends Clade> T getSecondClade(T clade1, T clade2) {
        if (clade1.getCladeInBits().nextSetBit(0) > clade2.getCladeInBits().nextSetBit(0)) {
            return clade1;
        } else {
            return clade2;
        }
    }

    public static int getFirstChild(Node vertex) {
        int minNr0 = getMinLeafNr(vertex.getChild(0));
        int minNr1 = getMinLeafNr(vertex.getChild(1));

        if (minNr0 < minNr1) {
            return 0;
        } else {
            return 1;
        }
    }

    public static int getSecondChild(Node vertex) {
        int minNr0 = getMinLeafNr(vertex.getChild(0));
        int minNr1 = getMinLeafNr(vertex.getChild(1));

        if (minNr0 > minNr1) {
            return 0;
        } else {
            return 1;
        }
    }

    private static int getMinLeafNr(Node vertex) {
        if (vertex.isLeaf())
            return vertex.getNr();
        else
            return vertex.getAllLeafNodes().stream().mapToInt(x -> x.getNr()).min().orElseThrow();
    }

    public static double getTimeSinceLastSplit(Node vertex) {
        double timeSinceLastSplit = Double.POSITIVE_INFINITY;

        for (Node child : vertex.getChildren()) {
            double currentTimeSinceLastSplit;
            if (child.isLeaf()) {
                currentTimeSinceLastSplit = vertex.getHeight() - child.getHeight();
            } else {
                currentTimeSinceLastSplit = getTimeSinceLastSplit(child);
            }
            timeSinceLastSplit = Math.min(currentTimeSinceLastSplit, timeSinceLastSplit);
        }

        return timeSinceLastSplit;
    }
}
