package ccd.model;

public class SBCCDClade extends Clade {
    public SBCCDClade(BitSet cladeInBits, AbstractCCD abstractCCD) {
        super(cladeInBits, abstractCCD);
    }

    @Override
    public CladePartition createCladePartition(Clade firstChildClade, Clade secondChildClade,
                                               boolean storeParent) {
        SBCCDClade[] partitioningClades = new SBCCDClade[]{(SBCCDClade) firstChildClade, (SBCCDClade) secondChildClade};

        CladePartition newPartition = new SBCCDCladePartition(this, partitioningClades);
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
