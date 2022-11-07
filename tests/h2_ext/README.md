# Extended heritability run

Test to run the heritability code for a longer period to see the full path to adaptation.
Previous test followed the first 100 generations of adaptation for a total 101 samples during the test period
We can't sample every generation for thousands of generations, so we'll sample every 20 gens for 2000 generations
for 100 samples. This should capture the adaptive walk: with 5% selection, it took ~1250 generations for adaptation to happen in my confirmation data.

I'm also not saving the relatedness matrix: it takes too much space, and I can recalculate it was a kinship matrix from the pedigree + phenotypes anyway