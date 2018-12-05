To align your images, create a subfolder within the find_roi directory. Number them as X ligation Y.tif, where each tif should have 4 channels but only 1 slice in Z, and X is constant and Y = 1...9999.

The ligations should be arranged in the order in which they are done, i.e.Primer N ligation 1, Primer N ligation 2, Primer N-1 ligation 1, Primer N-1 ligation 2, PRimer N-2 ligation 2 etc.

Then, run find_roi_stack as a script.

To generate barcodes from the aligned images of beads, run BeadSeq.m, giving it the base name to the ligation files. It will identify automatically how many ligations it has to process, but for debugging purposes you can set m in the for loop to go over a restricted number of values.

To generate barcodes from illumina sequencing data, run BeadBarcodeIdentification.m

To compare barcodes from Illumina and Solid, run CompareBarcodes.m