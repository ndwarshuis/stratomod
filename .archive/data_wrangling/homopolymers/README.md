# Perfect+imperfect homopolymer (generate separate annotations for A or T and G or C homopolymers, and make homopolymer length the 4th column)

 - Can build on JustinW’s script at: https://github.com/genome-in-a-bottle/genome-stratifications/blob/master/GRCh38/LowComplexity/GRCh38_SimpleRepeat_imperfecthomoplgt20_slop5.sh 
 - Run lines 11-19, except remove the selection of size > 20 (awk ‘$3-$2>20’) since we want to keep all sizes >3bp
 - Lines 21-34 take the union of A, T, C, and G homopolymer regions and add 5bp slop. We instead want to separately take the union of A+T and C+G homopolymers, and add the length to the 4th column before adding 5bp slop
 - See slack message at https://nist-oism.slack.com/archives/C02K2TVUPGW/p1637596400003600  for example merging for C+G homopolymers
 - Features most likely to be useful: C+G homopolymer length and A+T homopolymer length


## dependencies
- python 3.8.10
- bedtools => v2.27.1

``` bash
# findSimpleRegions_quad.py is from: https://opendata.nist.gov/pdrsrv/mds2-2190/GRCh38/LowComplexity/findSimpleRegions_quad.py
python findSimpleRegions_quad.py -p 20 -d 100000 -t 100000 -q 100000 GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GRCh38_SimpleRepeat_p20.bed

sed ‘s/^chr//’ GRCh38_SimpleRepeat_p20.bed | grep “^[0-9XY]” | grep -v ‘_’ | sed ‘s/^X/23/;s/^Y/24/’ | sort -k1,1n -k2,2n -k3,3n | sed ‘s/^23/X/;s/^24/Y/;s/^/chr/’ | bgzip -c > GRCh38_SimpleRepeat_homopolymer_gt20.bed.gz

python findSimpleRegions_quad.py -p 3 -d 100000 -t 100000 -q 100000 GCA_000001405.15_GRCh38_no_alt_analysis_set.fa GRCh38_SimpleRepeat_p3.bed

grep 'unit=C' GRCh38_SimpleRepeat_p3.bed | mergeBed -i stdin -d 1 > GRCh38_SimpleRepeat_imperfecthomopolgt3_C.bed

grep 'unit=G' GRCh38_SimpleRepeat_p3.bed | mergeBed -i stdin -d 1 > GRCh38_SimpleRepeat_imperfecthomopolgt3_G.bed

grep 'unit=A' GRCh38_SimpleRepeat_p3.bed | mergeBed -i stdin -d 1 > GRCh38_SimpleRepeat_imperfecthomopolgt3_A.bed

grep 'unit=T' GRCh38_SimpleRepeat_p3.bed | mergeBed -i stdin -d 1 > GRCh38_SimpleRepeat_imperfecthomopolgt3_T.bed

multiIntersectBed -i GRCh38_SimpleRepeat_imperfecthomopolgt20_C.bed \
	GRCh38_SimpleRepeat_imperfecthomopolgt20_G.bed \
	GRCh38_SimpleRepeat_imperfecthomopolgt20_A.bed \
	GRCh38_SimpleRepeat_imperfecthomopolgt20_T.bed |
    sed ‘s/^chr//’ |
    cut -f1-3 | grep “^[0-9XY]” | grep -v ‘_’ |
    sed ‘s/^/chr/’ |
	slopBed -i stdin -b 5 -g human.b38.genome |
    sed ‘s/^chr//’ |
	sed ‘s/^X/23/;s/^Y/24/’ |
	sort -k1,1n -k2,2n -k3,3n |
	sed ‘s/^23/X/;s/^24/Y/;s/^/chr/’ |
	mergeBed -i stdin |
	bgzip -c > GRCh38_SimpleRepeat_imperfecthomopolgt20_slop5.bed.gz

multiIntersectBed -i GRCh38_SimpleRepeat_imperfecthomopolgt3_C.bed GRCh38_SimpleRepeat_imperfecthomopolgt3_G.bed | sed 's/^chr//' | cut -f1-3 | grep "^[0-9XY]" | grep -v '_' | sed 's/^/chr/' | slopBed -i stdin -b 5 -g human.b38.genome | sed 's/^chr//' | sed 's/^X/23/;s/^Y/24/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^23/X/;s/^24/Y/;s/^/chr/' | mergeBed -i stdin | awk '{FS=OFS="\t"} {print $1,$2,$3,$3-$2-10}' > GRCh38_SimpleRepeat_imperfectCGhomopolgt3_slop5_size.bed


multiIntersectBed -i GRCh38_SimpleRepeat_imperfecthomopolgt3_A.bed GRCh38_SimpleRepeat_imperfecthomopolgt3_T.bed | sed 's/^chr//' | cut -f1-3 | grep "^[0-9XY]" | grep -v '_' | sed 's/^/chr/' | slopBed -i stdin -b 5 -g human.b38.genome | sed 's/^chr//' | sed 's/^X/23/;s/^Y/24/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^23/X/;s/^24/Y/;s/^/chr/' | mergeBed -i stdin | awk '{FS=OFS="\t"} {print $1,$2,$3,$3-$2-10}' > GRCh38_SimpleRepeat_imperfectAThomopolgt3_slop5_size.bed
```

