# EBDT
Python script for the EBDT algorithm

## Arguments
The following arguments allow parameterisation of the algorithm;

**-a** Kinase inhibition specificity threshold (default = 0.5)  
**-r** Phosphosite/kinase fingerprint ratio threshold (default = 0.5)  
**-p** Phosphosite being PDT of kinase probability threshold (default = 0.75)  
**-c** Specify the list of cell line file names separated by a comma. E.g.:  
       ```
       -c MCF7.xlsm,HL60.xlsm
       ```  
       Each specified cell line must appear in the listOfKinases.csv file. I.e. if "MCF7.xlsm" is specified, then there must be a row entry for "MCF7" under the cell line column.

## How to use
The algorithm is run by calling ebdt.py with at least the -c parameter detailing the list of cell line files to use, e.g.  
       ```
       python ebdt.py -c MCF7.xlsm,HL60.xlsm
       ```  
Additionally, user-defined thresholds can be applied to parameterise certain parts of the workflow. E.g. for specifying a kinase inhibition specificity threshold of 0.6:  
       ```
       python ebdt.py -c MCF7.xlsm,HL60.xlsm -a 0.6
       ```  

## Required Files
1. **Kinase inhibitor selectivity dataset file.** This must list the remaining kinase activity, for each kinase, after inhibitor treatment in vitro (normalised from 0 to 1, where 1 is no inhibition). This file must be included in a sub-directory *requiredData* and be named *kinaseInhibitionSpecificity.csv*. An example file is provided in the repository: *requiredData/kinaseInhibitionSpecificity.csv*.
2. **List of kinases for each cell line**
This file must list the the kinases known to be expressed in all the cell lines specified when calling the script.
3. **Phosphoproteomics dataset of kinase inhibitor effects.** This must be an Excel file consisting of two worksheets named *pvalue* and *fold* in which the p-Values and fold-change values are provided, respectively, for each phosphorylation site/compound pair. This file must be located in the same directory from which the script is being run. Column heading formats must strictly follow that of the provided example file *requiredData/MCF7.xlsm* - this is to ensure that compound names can be matched to those in the kinase inhibitor selectivity dataset file.
