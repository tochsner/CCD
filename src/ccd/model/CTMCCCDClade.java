package ccd.model;

public class CTMCCCDClade extends Clade {
    public CTMCCCDClade(BitSet cladeInBits, AbstractCCD abstractCCD) {
        super(cladeInBits, abstractCCD);
    }

    @Override
    public CladePartition createCladePartition(Clade firstChildClade, Clade secondChildClade,
                                               boolean storeParent) {
        CTMCCCDClade[] partitioningClades = new CTMCCCDClade[]{
                (CTMCCCDClade) Utils.getFirstClade(firstChildClade, secondChildClade),
                (CTMCCCDClade) Utils.getSecondClade(firstChildClade, secondChildClade)
        };


        CladePartition newPartition = new CTMCCCDCladePartition(this, partitioningClades);
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
