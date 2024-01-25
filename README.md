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
```

## Software usage
### Example command line:
```nifH_palmscan -search_pp protein.faa -all -nifH -model v16.nifH.pssm -ppout nifH.scanned.faa -fevout nifH.scanned.fev -report nifH.scanned.txt```
#### Type palmscan -help for option details.

## Reference
A. Babaian and R. C. Edgar (2022), Ribovirus classification by a polymerase barcode sequence, PeerJ. https://peerj.com/articles/14055/
R. C. Edgar et al. (2021), Petabase-scale sequence alignment catalyses viral discovery, Nature 2022 https://www.nature.com/articles/s41586-021-04332-2
