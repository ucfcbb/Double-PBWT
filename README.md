# Double-PBWT

Double-PBWT: Haplotype Switch-Point Detection by Double PBWT

# To Run

1. Change into the `/src/` directory.
2. Run the `make` command in a terminal from within the directory. 
`<your-home-path>/doublePBWT/src$ make`
This creates an exectuable `d2pbwt_overlap_genetic`.

Possible error on `makedepend` will require the user to install it(on Ubuntu) using the following command in a terminal:
`sudo apt-get install xutils-dev`

# Usage

Usage: `./d2pbwt_overlap_genetic <vcf-file> <genetic-mapping-file> <min-length-in-cM> <max-overlap-in-sites>`

The output file name should be specified by the user as shown in the example below:

## Example

`./d2pbwt_overlap_genetic ../example/vcf/test_chr22.vcf ../example/genmap/22_cMMap 2 4 > results`

The different columns of output file `results` are explained below in order:
`k2, k3, physical-site-k2, physical-site-k3, reference-individual-ID, switched-individual-ID, hap1-1, hap1-2, hap2-1, hap2-2, type`

* __k2__: Ending site(inclusive) of the first IBD match (found by lagging PBWT)
* __k3__: starting site of the leading match found by lagging PBWT. While the actual starting position of this IBD can be before this site, k3 is the switch point.
* __physical-site-k2__: Physical location corresponding to site k2
* __physical-site-k3__: Physical location corresponding to site k3
* __reference-individual-ID__: Individual ID for which the IBD is contiguous
* __switched-individual-ID__: Individual ID that has the switch
* __hap1-1, hap1-2__: The haplotypes that form the first IBD match. hap1-1 corresponds to reference-individual-ID, while hap1-2 corresponds to switched-individual-ID
* __hap2-1, hap2-2__: The haplotypes that form the second IBD match. hap2-1 corresponds to reference-individual-ID, while hap2-2 corresponds to switched-individual-ID
* __type__: __A__ represents alternating matches

# Citations
