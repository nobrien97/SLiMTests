# Get h2 from h2_hsfs

This job calculates heritability and variance components using sommer mmer functions, using the data from the h2_hsfs slim experiment. 

Some estimates didn't calculate due to no phenotypic variance. In future, maybe I should investigate values further from 1 to avoid populations getting stuck at 0 (due to accumulating massive alpha values)