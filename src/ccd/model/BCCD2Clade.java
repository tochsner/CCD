package ccd.model;

public class BCCD2Clade extends BCCDClade {
    public BCCD2Clade(BitSet cladeInBits, AbstractCCD abstractCCD) {
        super(cladeInBits, abstractCCD);
    }

    @Override
    public CladePartition createCladePartition(Clade firstChildClade, Clade secondChildClade,
                                               boolean storeParent) {
        BCCD2Clade[] partitioningClades = new BCCD2Clade[]{(BCCD2Clade) firstChildClade, (BCCD2Clade) secondChildClade};

        CladePartition newPartition = new BCCD2CladePartition(this, partitioningClades);
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
