
from cc3d import CompuCellSetup
        

from celldeath_v3Steppables import celldeath_v3Steppable

CompuCellSetup.register_steppable(steppable=celldeath_v3Steppable(frequency=1))



# from celldeath_v3Steppables import GrowthSteppable

# CompuCellSetup.register_steppable(steppable=GrowthSteppable(frequency=1))

CompuCellSetup.run()



# from cc3d import CompuCellSetup
        


# from celldeath_v2Steppables import ConstraintInitializerSteppable

# CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable(frequency=1))




# from celldeath_v2Steppables import GrowthSteppable

# CompuCellSetup.register_steppable(steppable=GrowthSteppable(frequency=1))


# CompuCellSetup.run()