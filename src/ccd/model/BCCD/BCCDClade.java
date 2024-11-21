package ccd.model.BCCD;

import ccd.model.AbstractCCD;
import ccd.model.BitSet;
import ccd.model.Clade;
import ccd.model.CladePartition;

public class BCCDClade extends Clade {
    public BCCDClade(BitSet cladeInBits, AbstractCCD abstractCCD) {
        super(cladeInBits, abstractCCD);
    }

    @Override
    public CladePartition createCladePartition(Clade firstChildClade, Clade secondChildClade,
                                               boolean storeParent) {
        BCCDClade[] partitioningClades = new BCCDClade[]{(BCCDClade) firstChildClade, (BCCDClade) secondChildClade};

        CladePartition newPartition = new BCCDCladePartition(this, partitioningClades);
        partitions.add(newPartition);

        childClades.add(partitioningClades[0]);
        childClades.add(partitioningClades[1]);
        if (storeParent) {
            partitioningClades[0].parentClades.add(this);
            partitioningClades[1].parentClades.add(this);
        }

        return newPartition;
    }
}
