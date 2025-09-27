# UCSC Annotation Data

* Human (*Homo sapiens*)  annotation BED files were downloaded via browser on Mon Mar 7 from the UCSC Browser with the following configurations:
    - `UCSC_canonical.bed.gz`
        - track: GENCODE V39
        - table: knownCanonical
        - region: genome
        - output format: selected fields from primary and related tables
        - file type returned: gzip compressed
        - Click "get output"
            - Select Fields from hg38.knownCanonical
                - chrom
                - chromStart
                - chromEnd
                - transcript
            - hg38.kgXref
                - geneSymbol
        - Click "get output"
    - `UCSC_exons.bed.gz`
        - track: GENCODE V39
        - table: knownGene
        - region: genome
        - output format: BED - browser extensible data
        - file type returned: gzip compressed
        - Click "get output"
            - Create one BED record per: Exons plus 0 bases at each end
    - `UCSC_introns.bed.gz`
        - track: GENCODE V39
        - table: knownGene
        - region: genome
        - output format: BED - browser extensible data
        - file type returned: gzip compressed
        - Click "get output"
            - Create one BED record per: Introns plus 0 bases at each end

## Modification of BED files

Modification of `UCSC_exons.bed.gz` and `UCSC_introns.bed.gz` to extract gene identifier and exon/intron info to create files `UCSC_exons_modif.bed` and `UCSC_introns.bed`:

    gzip -cd UCSC_exons.bed.gz | gawk '{split ($4,a,"_"); {print $1"\t"$2"\t"$3"\t"a[1]"\t"a[3]"\t"a[2]"\t"$6}}' > UCSC_exons_modif.bed

    gzip -cd UCSC_introns.bed.gz | gawk '{split ($4,a,"_"); {print $1"\t"$2"\t"$3"\t"a[1]"\t"a[3]"\t"a[2]"\t"$6}}' > UCSC_introns_modif.bed

Equality join of canonical and intron/exon modified tables by ENSEMBLID (Unix join command will not work unless both files are sorted by the column to join on) to create files `UCSC_exons_modif_canonical.bed` and `UCSC_introns_canonical.bed`.

    gzip -cd UCSC_canonical.bed.gz > UCSC_canonical.bed

    join -1 4 -2 4 <(sort -k4 UCSC_exons_modif.bed ) <(sort -k4 UCSC_canonical.bed) | gawk '{print $2"\t"$3"\t"$4"\t"$11"\t"$6"\t"$5}' | bedtools sort -i "-" > UCSC_exons_modif_canonical.bed

    join -1 4 -2 4 <(sort -k4 UCSC_introns_modif.bed ) <(sort -k4 UCSC_canonical.bed) | gawk '{print $2"\t"$3"\t"$4"\t"$11"\t"$6"\t"$5}' | bedtools sort -i "-" > UCSC_introns_modif_canonical.bed

    rm UCSC_canonical.bed

## CHECKSUMS (SHA-1)
- `UCSC_canonical.bed`: 69b967ffae859ec7fb9ffc02eb4800da0cc6f4a4
- `UCSC_exons.bed`: 409d2c6ca3dff0764ea3454e3db15f420142fa6e
- `UCSC_introns.bed`: cf276ae8c5714b1b09d612532a3879fa300c3daa
