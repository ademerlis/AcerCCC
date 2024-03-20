# How does an urbanized environment influence gene expression of *Acropora cervicornis*?

<img width="656" alt="Screen Shot 2023-08-04 at 11 47 45 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/db381d49-a962-4656-ad18-21efcd8a5c77">

From [Enochs et al. 2023](https://www.nature.com/articles/s41598-023-33467-7)


## Methods

**1. Field collections**

In June 2021, four genotypes of *A. cervicornis* were outplanted to the urban Coral City Camera site following collection from the in-situ coral nursery off Key Biscayne, FL. In late October/early November 2021, ~1 cm tissue samples were collected from each surviving outplant and an offshore in-situ nursery conspecific and immediately preserved in Zymo DNA/RNA Shield. 

![Acer color zipties at CCC June 2021](https://github.com/ademerlis/AcerCCC/assets/56000927/fa307ad4-e7ec-4225-b149-f7c8bd39edc5)
![Acer growth at CCC Nov 2021](https://github.com/ademerlis/AcerCCC/assets/56000927/828eb7dd-2f96-462b-a472-ef618a68dcc0)
Figure 2. Images from Dalton Hesley of initial outplants in June 2021 with colored zipties for genotype tracking, and different outplants after six months in the field (Nov 2021).

**2. Extractions and Sequencing**

Total RNA was extracted using the Zymo MagBead DNA/RNA extraction kit, then prepared as QuantSeq 3'mRNA-Seq cDNA libraries and sequenced on a NOVAseq S2 flow cell at the University of Miami Genomics Center.

**3. Bioinformatics**

Pipeline: FastQC -> TrimGalore (adapters + low-quality bp) -> TrimGalore (polyA tail) -> FastQC -> STAR -> DESeq2

See [this README](https://github.com/ademerlis/AcerCCC/blob/main/gene_expression/bioinformatics/README.md) for in-depth scripts of pipeline.

