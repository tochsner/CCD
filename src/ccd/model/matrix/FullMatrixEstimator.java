package ccd.model.matrix;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import ccd.model.HeightSettingStrategy;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class FullMatrixEstimator extends MatrixEstimator {
    @Override
    public TreeMatrixConfiguration estimateConfiguration(MatrixCCD ccd) {
        List<Pair<Integer, Integer>> allDistances = new ArrayList<>();

        for (int i = 0; i < ccd.getNumberOfLeaves(); i++) {
            for (int j = i + 1; j < ccd.getNumberOfLeaves(); j++) {
                allDistances.add(new Pair<>(i, j));
            }
        }

        TreeMatrixConfiguration configuration = CubeUtils.getConfiguration(
                allDistances,
                ccd.getMAPTree()
        );
        configuration.allTreesAreCompatible = true;
        return configuration;
    }
}
