# nifH_palmscan

## nifH_Palmscan algorithm
nifH_palmscan is software designed to scan the palmprint of the nitrogenase iron protein (a segment with well-conserved catalytic motifs) from DNA/protein file. 
nifH-Palmscan is developed based on the framework of Palmscan (A. Babaian and R. C. Edgar, 2022).

## Repository layout
```
palmscan/
  src/               # Source code (C++)
  v16.nifH.pssm      # PSSM matrices
  nifH_palmscan      # Binary file
  protein.faa        # Demo for test
```

## OS Requirements
```
The developmental version of the package has been tested on the following architecture and operating systems:
Linux: Red Hat 4.8.5-44
Architecture: x86_64
```

## Installation Guide
```
yum install ccache gcc-c++ make git glibc-static -y
git clone https://github.com/JaphethLin/nifH_palmscan.git
cd palmscan/src
make
```

## Software usage
### Example command line:
```nifH_palmscan -search_pp protein.faa -all -nifH -model v16.nifH.pssm -ppout nifH.scanned.faa -fevout nifH.scanned.fev -report nifH.scanned.txt```
#### Type palmscan -help for option details.

## Results
```
nifH.scanned.faa      # All scanned nifH palmprint sequence fragments
nifH.scanned.fev      # Details of scoring
nifH.scanned.txt      # Sequences aligned to the model
```

## Reference

A. Babaian and R. C. Edgar, Ribovirus classification by a polymerase barcode sequence, PeerJ (2022) https://peerj.com/articles/14055

R. C. Edgar et al., Petabase-scale sequence alignment catalyses viral discovery, Nature (2021) https://www.nature.com/articles/s41586-021-04332-2
