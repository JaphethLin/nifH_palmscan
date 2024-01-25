# nifH_palmscan

## nifH_Palmscan algorithm
nifH-Palmscan developed based on the framework of Palmscan (A. Babaian and R. C. Edgar, 2022) 

nifH_palmscan modeling the palmprint of the nitrogenase iron protein — a segment of the sub-domain with well-conserved catalytic motifs — for scanning it from DNA/protein file.

## Repository layout
```
palmscan/
  src/               # Source code (C++)
  v16.nifH.pssm      # PSSM matrices
  nifH_palmscan      # Binary file
```

## Software usage
### Example command line:
```nifH_palmscan -search_pp protein.faa -all -nifH -model v16.nifH.pssm -ppout nifH.scanned.faa -fevout nifH.scanned.fev -report nifH.scanned.txt```
### Type palmscan -help for option details.

## Reference
A. Babaian and R. C. Edgar (2022), Ribovirus classification by a polymerase barcode sequence, PeerJ. https://peerj.com/articles/14055/
R. C. Edgar et al. (2021), Petabase-scale sequence alignment catalyses viral discovery, Nature 2022 https://www.nature.com/articles/s41586-021-04332-2
