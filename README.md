# VME: Visual Microbiome Explorer - A Tool for Visual Exploration of Microbial Communities

VME is a highly interactive tool for visual analysis of microbial communities. VME employ complete linked views and full user interaction between the tool and the data. Moreover, integrates metadata in the exploration process, speeding up the analysis, and gaining insights from the data.

## Requirements:
 * R >= 4.0.3
 * RStudio >= 1.3 (optional)
 * installapps.R (R packages used in VME)
 * Docker (recommended). Optional for running Docker image of VME.

## Folders:
 * app
   * app.R shiny code
   * processingfunctions.R functions code

 * datasets
   * Endesfelder
     * otu_data.tsv
     * metadata.txt
     * 97_otus.tree
     * communities1.txt
   * Wegner
     * data_set.tsv
     * metadata.txt
   * Bazanella
     * bazanella_mth12.tsv
     * Metafile_Month12.txt
     * Rooted_tree_Month12.nwk

## Installation
**Install shiny and the R packages to run VME**
clone folder https://github.com/charlos1204/VME.git and run script installapps.R 
```bash
git clone https://github.com/charlos1204/VME.git
cd VME
Rscript installapps.R
```
or run R and

```r
source("installapps.R")
```

**Run VME app**
Edit app.R and change line 25 to VME path folder.
Change to VME folder and run R.

```r
library(shiny)
shiny::runApp('app')
```
In RStudio open the app.R file and click in Run app button.

**Run with Docker (recommended)**
Pull the Docker image:
```bash
docker pull charlos1204/vme:latest
```

Run Docker image:
```bash
docker run -p 3838:3838 charlos1204/vme
```
Copy link http://0.0.0.0:3838 into a browser
