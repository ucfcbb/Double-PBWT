# Double-PBWT

Double-PBWT: Haplotype Switch-Point Detection by Double PBWT

# To Run

1) Change into the /src/ directory.
2) Run the make command in the directory from a terminal.
<your-home-path>/doublePBWT/src$ make

This creates an exectuable d2pbwt_overlap_genetic.

# Usage

Usage: ./d2pbwt_overlap_genetic <vcf-file> <genetic-mapping-file> <min-length-in-cM> <max-overlap-in-sites>

The output file name should be specified by the user as shown in the example below:

## Example

./d2pbwt_overlap_genetic ../example/vcf/test_chr22.vcf ../example/genmap/22_cMMap 2 4 > results
