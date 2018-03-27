This folder contains the data used for RNA-seq experiments simulations. 

The `iso_info-tab` file is a matrix containing isoform information. The specified columns are:
  
  - `transcript_id`: Transcript/isoform identificator.
  - `gene_id`: Gene identificator
  - `numb_iso`: Amount of isoforms of the gene.
  - `group`: Gene group definid according to the number of isoforms per gene. Group ''0', one isoform per gene; group 'I', two to four isoforms; group 'II', five to nine isoforms; group 'III', more than nine isoforms.
  - `meanC1`: Mean gene expression observed in the real dataset for the Control condition.
  - `meanC2`: Mean gene expression observed in the real dataset for the Treatment condition.
  - `mean`: Mean gene expression.
  - `varC1`: Averaged observed gene variance in the real dataset for the Control condition.
  - `varC2`: Averaged observed gene variance in the real dataset for the Treatment condition.
  - `var`: Averaged observed gene variance.
  - `iso_meanC1`: Mean isoform expression observed in the real dataset for the Control condition.
  - `iso_meanC2`: Mean isoform expression observed in the real dataset for the Treatment condition.
  - `iso_mean`: Mean isoform expression.
  - `iso_varC1`: Averaged observed isoform variance in the real dataset for the Control condition.
  - `iso_varC2`: Averaged observed isoform variance in the real dataset for the Treatment condition.
  - `iso_var`: Averaged observed isoform variance.  
  - `ratiosC1`: Mean isoform proportion (relative expression) observed in the real dataset for the Control condition.
  - `ratiosC2`: Mean isoform proportion observed in the real dataset for the Treatment condition.
  - `DS`: Logical indicating if the gene is simulated as 'DS'.
  - `DIEDS`: Logical indicating if the gene is simulated as 'DIEDS'.
  - `DSgroup`: If the gene is considered in the 'DS' simulation group, this column indicates the simulated change. Else, contains zero values.
  - `DIEDSgroup`: If the gene is considered in the 'DIEDS' simulation group, this column indicates the simulated change. Else, contains zero values.
  - `meanC1Sim`: Simulated mean isoform expression for the Control condition.
  - `meanC2Sim`: Simulated mean isoform expression for the Treatment condition.
  - `var1Sim`: Simulated isoform expression variance for the Control condition.
  - `varC2Sim`: Simulated isoform expression variance for the Treatment condition.
  - `mayor`: Logical indicating if the isoform is the major (most expressed) isoform of the gene.
  - `updown`: In DS and DIEDS genes indicates the direction of the change.
  
The `gene2transc.txt` file contains the gene-isoform relationship used here.

Each folder called `ScenarioY` (Y={1,2,3}) contains ten folders called`simX` (X=1,2,..,10) which contain the simulation files obtained for each experiment replication. The `iso_cm_sim.RData` R dataset contains the simualted expression matrix at the isoform level. The `X_sim_iso.results` are the RSEM quantification files with the simulated profiles for each sample.
  
 
