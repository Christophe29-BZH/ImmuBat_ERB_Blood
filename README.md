# Blood Leukocytes from the Egyptian Rousette Bat
This repository contains scripts for the single cell analysis of the Egyptian Rousette Bat (ERB) blood leukocytes and generating figures used in the Cell Report paper.


Blood samples were obtained from a total of 9 bats across different ages. The analysis aims at identifying the different cell types / states present in ERB peripheral blood and their abundancy variations upon ageing.

## Data Accessibility
Sequencing data and count matrices are stored on [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183925).

Gene expression data can be directly browsed online on the [Infection Atlas](https://infection-atlas.org/Immubat-ERB-Blood/).

No manual download of data is required before running the scripts. 
The first script will automatically downloads the count matrices from GEO.

## Repository Structure
**data:** Directory containing additional raw data. This includes the output of Souporcell to demultiplex the libraries, flow cytometry data obtained when processing the samples (*Data_FlowCytometry_FigS2.csv*) and general information about the composition of the libraries (*Library_info.csv*). Data downloaded from GEO will be stored in this directory too.

**analysis:** Directory containing the R and python scripts.

**outputs:** This directory is created when running the scripts. It will contain the processed data and different results including differentially expressed gene tables, annotated Seurat objects, plots...

## Code Execution
#### 1- Download repository
Option 1: Download manually the repository as a ZIP archive and extract it locally on your computer

Option 2: Clone the repository
```shell
git clone https://github.com/Christophe29-BZH/ImmuBat_ERB_Blood.git
cd ImmuBat_ERB_Blood/analysis
```


#### 2- Install R and python dependencies 
See Dependencies section.


#### 3- Run sequentially the scripts in the **analysis** directory 
Make sure to set the **analysis** directory as the working directory when running the scripts.

## Dependencies
List of R and python packages necessary to run the scripts.

#### Required libraries for R
- R                 4.0.3
- Seurat            4.0.0
- stringr           1.4.0
- dplyr             1.0.7
- ggplot2           3.3.5
- cowplot           1.1.1
- igraph            1.2.6
- leidenbase        0.1.3
- reshape2          1.4.4 
- grid              4.0.3
- tidyr             1.1.4
- RColorBrewer      1.1-2
- data.table        1.14.0
- stringi           1.7.5
- ggpubr            0.4.0

#### Required libraries for python
- python            3.9.7
- sccoda            0.1.7
- notebook          6.4.8
