# scHaplotyper: Haplotype construction and visualization for genetic diagnosis using single cell DNA sequencing data

![](https://github.com/yzqheart/FiguresForRef/blob/master/For.scHaplotyper.readme.ref.fig1.jpg)<br>

scHaplotyper performs haplotyping for genetic diagnosis using single cell DNA sequencing data. First, trio father-mother-child and father-mother-embryo “core families” are constructed, and the affected child and embryos are phased to paternal and maternal inheritance. The affected child-inherited paternal and maternal haplotypes are shown as red and yellow, respectively. Second, the haplotypes of the embryos are reconstructed by regarding the affected child as a reference. The horizontal line indicates the mutation position on the chromosome. The embryo is diagnosed as carrying paternal and/or maternal mutations when the embryo inherits the same paternal and/or maternal haplotype block as the affected child at the mutation position on the chromosome.<br>

## Platform and prerequisite
* Linux<br>
* bcftools<br>
* Perl [SVG](https://metacpan.org/release/SVG) module is required.<br>
To install perl SVG, simply copy and paste either of the commands in to your terminal:<br>
	* cpanm:<br>
		`cpanm SVG`<br>
    
	* CPAN shell:<br>
		`perl -MCPAN -e shell`<br>
		`install SVG`<br>
Test SVG install: `perldoc SVG`<br>    
## Obtain and install
Download the latest source release on github.<br>

Unpack the file and change directory to the package directory as follows:<br>

`unzip scHaplotyper-master.zip`<br>

`cd scHaplotyper-master`<br>

Configure the package:<br>

`chmod +x scHaplotyper bin/*`<br>

## Data requirements
scHaplotyper requires:<br>
* sample information file indicates identity of each sample and known disease allea carrier status of the diagnosis family.<br>
* mutation information contains mutation type and position of the gene.<br>
* joint called vcf file of the family.<br>

See the above [example](https://github.com/yzqheart/scHaplotyper/tree/master/example) to explore the example inputs provided.

## Usage and example
### Usage:<br>
`./scHaplotyper SampleInfoFile MutationInfoFile VcfFile`<br>

### Example:<br>
`./scHaplotyper example/Case2.SampleInfo.csv example/Case2.MutationInfo.csv example/Case2.vcf`<br>

