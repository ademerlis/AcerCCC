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

See [this README](https://github.com/ademerlis/AcerCCC/blob/main/bioinformatics/README.md) for in-depth scripts of pipeline.

## Results

**Sequencing Quality**

### **1. Raw reads**

<img width="840" alt="Screen Shot 2023-08-21 at 9 45 47 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/9c737c7a-aaef-453c-8b4e-cb1a310a4f95">

<img width="849" alt="Screen Shot 2023-08-21 at 9 46 03 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/9a3d3c1b-d545-4a50-a8e9-34a773d67f62">

<img width="861" alt="Screen Shot 2023-08-21 at 9 46 16 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/b78f7ae6-a5f3-4536-a209-d90b3d814e4d">

<img width="856" alt="Screen Shot 2023-08-21 at 9 46 33 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/286b67ae-f321-4409-a84f-29d261dbd452">

<img width="864" alt="Screen Shot 2023-08-21 at 9 46 52 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/9fe2db6f-e019-4d09-b0bc-64eda2c2c94d">

<img width="869" alt="Screen Shot 2023-08-21 at 9 47 08 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/9db6bb5a-a265-44f1-9f27-8df440bf046f">

_______________

### **2. Trimmed reads (adapters script)**

<img width="856" alt="Screen Shot 2023-08-21 at 9 50 51 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/da695c60-032c-40f7-ab5a-1d7523353f6b">

<img width="863" alt="Screen Shot 2023-08-21 at 10 07 17 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/48beae46-961b-4f3e-87f2-2dd9350c76f9">

<img width="848" alt="Screen Shot 2023-08-21 at 10 03 35 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/6e6a3e89-29eb-478c-9e8c-9709bd1adcb9">

<img width="848" alt="Screen Shot 2023-08-21 at 10 03 49 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/e4fe4e38-6244-4e6d-9473-66d77a244777">

<img width="861" alt="Screen Shot 2023-08-21 at 10 03 58 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/9499f4c6-a330-4ed6-b845-c0903b8dc155">

_______________

### **3. Trimmed reads (polyA tail script)**

<img width="847" alt="Screen Shot 2023-08-21 at 9 51 25 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/10ff8d49-ffbc-49bd-8cbb-7a47dbea65ee">

<img width="850" alt="Screen Shot 2023-08-21 at 9 52 15 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/ca2400e3-380b-4621-9c1e-2cc4ce7756a3">

<img width="849" alt="Screen Shot 2023-08-21 at 10 07 39 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/e236d969-8bf1-47d1-ba72-03a6be766bbd">

<img width="864" alt="Screen Shot 2023-08-21 at 10 07 51 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/76468320-60b1-4d16-8817-ff716d8afdcb">

<img width="855" alt="Screen Shot 2023-08-21 at 10 08 00 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/30957b3a-42e2-47fe-a694-5c0076a0f016">

<img width="857" alt="Screen Shot 2023-08-21 at 10 08 10 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/01f0bfce-37ce-4160-825f-aee3e640e806">

<img width="856" alt="Screen Shot 2023-08-21 at 10 08 22 AM" src="https://github.com/ademerlis/AcerCCC/assets/56000927/cc8a176b-283f-49af-8cea-77bbe8e7af1e">

