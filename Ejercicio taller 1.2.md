# Tarea Sesión 
Joaquín Ignacio Torres Núñez

**_Bioinformatica y analisis genómicos reproducibles_**

## **Ejercicio 1: Ejercicio: abre el el editor de Markdown de tu preferencia y escribe un texto en formato Markdown de manera que quede igual que los tres primeros puntos de Preparing the environment, cleaning the data for Stacks (incluyendo ese subtítulo). No es necesario poner los colores, pero si quieres, cool.**

# 4. Running the pipeline
## 4.1 Clean the data
In a typical analysis, data will be received from an Illumina sequencer, or some other type of sequencer as FASTQ files. The first requirement is to demultiplex, or sort, the raw data to recover the individual samples in the Illumina library. While doing this, we will use the Phred scores provided in the FASTQ files to discard sequencing reads of low quality. These tasks are accomplished using the <font color="green">process_radtags</font> program. 
![alt text](https://catchenlab.life.illinois.edu/stacks/manual/process_radtags.png)
Some things to consider when running this program: 
- <font color="green">process_radtags</font> can handle both single-end or paired-end Illumina sequencing.
- The raw data can be compressed, or gzipped (files end with a "```.gz```" suffix).
- You can supply a list of barcodes, or indexes, to <font color="green">process_radtags</font>  in order for it to demultiplex your samples. These barcodes can be single-end barcodes or combinatorial barcodes (pairs of barcodes, one on each of the paired reads). Barcodes are specified, one per line (or in tab separated pairs per line), in a text file.
    - If, in addition to your barcodes, you also supply a sample name in an extra column within the barcodes file, <font color="green">process_radtags</font>  will name your output files according to sample name instead of barcode.
    - If you supply the same sample name for multiple barcodes, reads containing those barcodes will be consolidated into a single output file, merging them.
- If you believe your reads may contain adapter contamination <font color="green">process_radtags</font>  can filter it out.
- You can supply the restriction enzyme used to construct the library. In the case of double-digest RAD, you can supply both restriction enzymes.
- If instructed, (-r command line option), <font color="green">process_radtags</font> will correct barcodes and restriction enzyme sites that are within a certain distance from the true barcode or restriction enzyme cutsite.
- By default, <font color="green">process_radtags</font>  will identify and filter reads containing runs of poly-Gs. These runs are often indicative of synethesis termination in two-channel sequencing platforms &mdash they represent the absence of sequence once the DNA molecule "runs out" during sequencing. The presence of these poly-G runs is platform dependent, thus <font color="green">process_radtags</font> will use the FASTQ headers to identify the correct platorm whenever possible and adjust filters accordingly. This behavior can be modified in the advanced options. 
### 4.1.1 Understanding barcodes/indexes and specifying the barcode type
Genotype by sequencing libraries sample the genome by selecting DNA adjacent to one or more restriction enzyme cutsites. By reducing the amount of total DNA sampled, most researchers will multiplex many samples into one molecular library. Individual samples are demarcated in the library by ligating an oligo barcode onto the restriction enzyme-associated DNA for each sample. Alternatively, an index barcode is used, where the barcode is located upstream of the sample DNA within the sequencing adaptor. Regardless of the type of barcode used, after sequencing, the data must be demultiplexed so the samples are again separated. The process_radtags program will perform this task, but we much specify the type of barcodes used, and where to find them in the sequencing data.
There are a number of different configurations possible, each of them is detailed below. 
1. If your data are **single-end** or **paired-end**, with an inline barcode present only on the single-end (marked in red):
```
@HWI-ST0747:188:C09HWACXX:1:1101:2968:2083 1:N:0:
TTATGATGCAGGACCAGGATGACGTCAGCACAGTGCGGGTCCTCCATGGATGCTCCTCGGTCGTGGTTGGGGGAGGAGGCA
+
@@@DDDDDBHHFBF@CCAGEHHHBFGIIFGIIGIEDBBGFHCGIIGAEEEDCC;A?;;5,:@A?=B5559999B@BBBBBA
@HWI-ST0747:188:C09HWACXX:1:1101:2863:2096 1:N:0:
TTATGATGCAGGCAAATAGAGTTGGATTTTGTGTCAGTAGGCGGTTAATCCCATACAATTTTACACTTTATTCAAGGTGGA
+
CCCFFFFFHHHHHJJGHIGGAHHIIGGIIJDHIGCEGHIFIJIH7DGIIIAHIJGEDHIDEHJJHFEEECEFEFFDECDDD
@HWI-ST0747:188:C09HWACXX:1:1101:2837:2098 1:N:0:
GTGCCTTGCAGGCAATTAAGTTAGCCGAGATTAAGCGAAGGTTGAAAATGTCGGATGGAGTCCGGCAGCAGCGAATGTAAA
```
Then you can specify the <code>--inline_null</code> flag to <font color="green">process_radtags</font> . This is also the default behavior and the flag can be ommitted in this case. 

2. If your data are **single-end** or **paired-end**, with a single index barcode (in blue):
```
@9432NS1:54:C1K8JACXX:8:1101:6912:1869 1:N:0:ATGACT
TCAGGCATGCTTTCGACTATTATTGCATCAATGTTCTTTGCGTAATCAGCTACAATATCAGGTAATATCAGGCGCA
+
CCCFFFFFHHHHHJJJJJJJJIJJJJJJJJJJJHIIJJJJJJIJJJJJJJJJJJJJJJJJJJGIJJJJJJJHHHFF
@9432NS1:54:C1K8JACXX:8:1101:6822:1873 1:N:0:ATGACT
CAGCGCATGAGCTAATGTATGTTTTACATTCCAGAAAGAGAGCTACTGCTGCAGGTTGTGATAAAATAAAGTAAGA
+
B@@FFFFFHFFHHJJJJFHIJHGGGHIJIIJIJCHJIIGGIIIGGIJEHIJJHII?FFHICHFFGGHIIGG@DEHH
@9432NS1:54:C1K8JACXX:8:1101:6793:1916 1:N:0:ATGACT
TTTCGCATGCCCTATCCTTTTATCACTCTGTCATTCAGTGTGGCAGCGGCCATAGTGTATGGCGTACTAAGCGAAA
+
@C@DFFFFHGHHHGIGHHJJJJJJJGIJIJJIGIJJJJHIGGGHGII@GEHIGGHDHEHIHD6?493;AAA?;=;=
```

Then you can specify the <code>--index_null</code> flag to <font color="green">process_radtags</font> . 

3. If your data are **single-end** with both an inline barcode (in red) and an index barcode (in blue):
```
@9432NS1:54:C1K8JACXX:8:1101:6912:1869 1:N:0:ATCACG
TCACGCATGCTTTCGACTATTATTGCATCAATGTTCTTTGCGTAATCAGCTACAATATCAGGTAATATCAGGCGCA
+
CCCFFFFFHHHHHJJJJJJJJIJJJJJJJJJJJHIIJJJJJJIJJJJJJJJJJJJJJJJJJJGIJJJJJJJHHHFF
@9432NS1:54:C1K8JACXX:8:1101:6822:1873 1:N:0:ATCACG
GTCCGCATGAGCTAATGTATGTTTTACATTCCAGAAAGAGAGCTACTGCTGCAGGTTGTGATAAAATAAAGTAAGA
+
B@@FFFFFHFFHHJJJJFHIJHGGGHIJIIJIJCHJIIGGIIIGGIJEHIJJHII?FFHICHFFGGHIIGG@DEHH
@9432NS1:54:C1K8JACXX:8:1101:6793:1916 1:N:0:ATCACG
GTCCGCATGCCCTATCCTTTTATCACTCTGTCATTCAGTGTGGCAGCGGCCATAGTGTATGGCGTACTAAGCGAAA
+
@C@DFFFFHGHHHGIGHHJJJJJJJGIJIJJIGIJJJJHIGGGHGII@GEHIGGHDHEHIHD6?493;AAA?;=;=
```
Then you can specify the <code>--inline_index</code> flag to <font color="green">process_radtags</font> .

4. If your data are **paired-end** with an inline barcode on the single-end (in red) and an index barcode (in blue):
```
@9432NS1:54:C1K8JACXX:7:1101:5584:1725 1:N:0:CGATGT
ACTGGCATGATGATCATAGTATAACGTGGGATACATATGCCTAAGGCTAAAGATGCCTTGAAGCTTGGCTTATGTT
+
#1=DDDFFHFHFHIFGIEHIEHGIIHFFHICGGGIIIIIIIIAEIGIGHAHIEGHHIHIIGFFFGGIIIGIIIEE7
@9432NS1:54:C1K8JACXX:7:1101:5708:1737 1:N:0:CGATGT
TTCGACATGTGTTTACAACGCGAACGGACAAAGCATTGAAAATCCTTGTTTTGGTTTCGTTACTCTCTCCTAGCAT
+
#1=DFFFFHHHHHJJJJJJJJJJJJJJJJJIIJIJJJJJJJJJJIIJJHHHHHFEFEEDDDDDDDDDDDDDDDDD@
```
```
@9432NS1:54:C1K8JACXX:7:1101:5584:1725 2:N:0:CGATGT
AATTTACTTTGATAGAAGAACAACATAAGCCAAGCTTCAAGGCATCTTTAGCCTTAGGCATATGTATCCCACGTTA
+
@@@DFFFFHGHDHIIJJJGGIIIEJJJCHIIIGIJGGEGGIIGGGIJIJIHIIJJJJIJJJIIIGGIIJJJIICEH
@9432NS1:54:C1K8JACXX:7:1101:5708:1737 2:N:0:CGATGT
AGTCTTGTGAAAAACGAAATCTTCCAAAATGCTAGGAGAGAGTAACGAAACCAAAACAAGGATTTTCAATGCTTTG
+
C@CFFFFFHHHHHJJJJJJIJJJJJJJJJJJJJJIJJJHIJJFHIIJJJJIIJJJJJJJJJHGHHHHFFFFFFFED
```
 Then you can specify the <code>--inline_index</code> flag to <font color="green">process_radtags</font> .

5. If your data are **paired-end** with indexed barcodes on the single and paired-ends (in blue):
```
@9432NS1:54:C1K8JACXX:7:1101:5584:1725 1:N:0:ATCACG+CGATGT
ACTGGCATGATGATCATAGTATAACGTGGGATACATATGCCTAAGGCTAAAGATGCCTTGAAGCTTGGCTTATGTT
+
#1=DDDFFHFHFHIFGIEHIEHGIIHFFHICGGGIIIIIIIIAEIGIGHAHIEGHHIHIIGFFFGGIIIGIIIEE7
@9432NS1:54:C1K8JACXX:7:1101:5708:1737 1:N:0:ATCACG+CGATGT
TTCGACATGTGTTTACAACGCGAACGGACAAAGCATTGAAAATCCTTGTTTTGGTTTCGTTACTCTCTCCTAGCAT
+
#1=DFFFFHHHHHJJJJJJJJJJJJJJJJJIIJIJJJJJJJJJJIIJJHHHHHFEFEEDDDDDDDDDDDDDDDDD@
```
```
@9432NS1:54:C1K8JACXX:7:1101:5584:1725 2:N:0:ATCACG+CGATGT
AATTTACTTTGATAGAAGAACAACATAAGCCAAGCTTCAAGGCATCTTTAGCCTTAGGCATATGTATCCCACGTTA
+
@@@DFFFFHGHDHIIJJJGGIIIEJJJCHIIIGIJGGEGGIIGGGIJIJIHIIJJJJIJJJIIIGGIIJJJIICEH
@9432NS1:54:C1K8JACXX:7:1101:5708:1737 2:N:0:ATCACG+CGATGT
AGTCTTGTGAAAAACGAAATCTTCCAAAATGCTAGGAGAGAGTAACGAAACCAAAACAAGGATTTTCAATGCTTTG
+
C@CFFFFFHHHHHJJJJJJIJJJJJJJJJJJJJJIJJJHIJJFHIIJJJJIIJJJJJJJJJHGHHHHFFFFFFFED
```
Then you can specify the <code>--index_index</code> flag to <font color="green">process_radtags</font> .

6. If your data are **paired-end** with inline barcodes on the single and paired-ends (in red):
```
@9432NS1:54:C1K8JACXX:7:1101:5584:1725 1:N:0:
ACTGGCATGATGATCATAGTATAACGTGGGATACATATGCCTAAGGCTAAAGATGCCTTGAAGCTTGGCTTATGTT
+
#1=DDDFFHFHFHIFGIEHIEHGIIHFFHICGGGIIIIIIIIAEIGIGHAHIEGHHIHIIGFFFGGIIIGIIIEE7
@9432NS1:54:C1K8JACXX:7:1101:5708:1737 1:N:0:
TTCGACATGTGTTTACAACGCGAACGGACAAAGCATTGAAAATCCTTGTTTTGGTTTCGTTACTCTCTCCTAGCAT
+
#1=DFFFFHHHHHJJJJJJJJJJJJJJJJJIIJIJJJJJJJJJJIIJJHHHHHFEFEEDDDDDDDDDDDDDDDDD@```
```@9432NS1:54:C1K8JACXX:7:1101:5584:1725 2:N:0:
AATTTACTTTGATAGAAGAACAACATAAGCCAAGCTTCAAGGCATCTTTAGCCTTAGGCATATGTATCCCACGTTA
+
@@@DFFFFHGHDHIIJJJGGIIIEJJJCHIIIGIJGGEGGIIGGGIJIJIHIIJJJJIJJJIIIGGIIJJJIICEH
@9432NS1:54:C1K8JACXX:7:1101:5708:1737 2:N:0:
AGTCTTGTGAAAAACGAAATCTTCCAAAATGCTAGGAGAGAGTAACGAAACCAAAACAAGGATTTTCAATGCTTTG
+
C@CFFFFFHHHHHJJJJJJIJJJJJJJJJJJJJJIJJJHIJJFHIIJJJJIIJJJJJJJJJHGHHHHFFFFFFFED
```
Then you can specify the <code>--inline_inline</code> flag to <font color="green">process_radtags</font> .   

### 4.1.2 Specifying the barcodes
The <font color="green">process_radtags</font>  program will demultiplex data if it is told which barcodes/indexes to expect in the data. It will also properly name the output files if the user specifies how to translate a particular barcode to a specific output file neam. This is done with the barcodes file, which we provide to <font color="green">process_radtags</font>  program. The barcode file is a very simple format — one barcode per line; if you want to rename the output files, the sample name prefix is provided in the second column. 
```
cat barcodes_lane3
CGATA<tab>sample_01
CGGCG     sample_02
GAAGC     sample_03
GAGAT     sample_04
TAATG     sample_05
TAGCA     sample_06
AAGGG     sample_07
ACACG     sample_08
ACGTA     sample_09
```
The sample names can be whatever is meaningful for your project: 
```
more barcodes_run01_lane01
CGATA<tab>spruce_site_12-01
CGGCG     spruce_site_12-02
GAAGC     spruce_site_12-03
GAGAT     spruce_site_12-04
TAATG     spruce_site_06-01
TAGCA     spruce_site_06-02
AAGGG     spruce_site_06-03
ACACG     spruce_site_06-04
```
Combinatorial barcodes are specified, one per column, separated by a tab: 
```
cat barcodes_lane07
CGATA<tab>ACGTA<tab>sample_01
CGGCG     ACGTA     sample_02
GAAGC     ACGTA     sample_03
GAGAT     ACGTA     sample_04
CGATA     TAGCA     sample_05
CGGCG     TAGCA     sample_06
GAAGC     TAGCA     sample_07
GAGAT     TAGCA     sample_08
```
Merging the first three sets of combinatorial barcodes into a single output file:
```
cat barcodes_lane07
CGATA<tab>ACGTA<tab>sample_01
CGGCG     ACGTA     sample_01
GAAGC     ACGTA     sample_01
GAGAT     ACGTA     sample_02
CGATA     TAGCA     sample_03
```

1.If you don’t want <font color="green">process_radtags</font>  to rename your samples, simply do not specify the last column in the barcodes file, and the output files will instead be named after the barcode.

2.Often, sequencing centers will return data from indexed libraries already demultiplexed. In this case, omit the barcodes file and <code>process_radtags</code> will not attempt to demulitplex the data, but can still be used to clean the data.

### 4.1.3 Running process_radtags
Here is how single-end data received from an Illumina sequencer might look: 
```
ls ./raw
lane3_NoIndex_L003_R1_001.fastq.gz  lane3_NoIndex_L003_R1_006.fastq.gz  lane3_NoIndex_L003_R1_011.fastq.gz
lane3_NoIndex_L003_R1_002.fastq.gz  lane3_NoIndex_L003_R1_007.fastq.gz  lane3_NoIndex_L003_R1_012.fastq.gz
lane3_NoIndex_L003_R1_003.fastq.gz  lane3_NoIndex_L003_R1_008.fastq.gz  lane3_NoIndex_L003_R1_013.fastq.gz
lane3_NoIndex_L003_R1_004.fastq.gz  lane3_NoIndex_L003_R1_009.fastq.gz
lane3_NoIndex_L003_R1_005.fastq.gz  lane3_NoIndex_L003_R1_010.fastq.gz
```
Then you can run <font color="green">process_radtags</font>  in the following way:

```
Process_radtags -p ./raw/ -o ./samples/ -b ./barcodes/barcodes_lane3 \
                  -e sbfI -r -c -q
```
I specify the directory containing the input files, <code>./raw</code>, the directory I want <font color="green">process_radtags</font>  to enter the output files, <code>./samples</code>, and a file containing my barcodes, <code>./barcodes/barcodes_lane3</code>, along with the restrction enzyme I used and instructions to clean the data and correct barcodes and restriction enzyme cutsites (```-r```, ```-c```, ```-q```).

Here is a more complex example, using paired-end double-digested data (two restriction enzymes) with combinatorial barcodes, and gzipped input files. Here is what the raw Illumina files may look like:  
```
ls ./raw
GfddRAD1_005_ATCACG_L007_R1_001.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_001.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_002.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_002.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_003.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_003.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_004.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_004.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_005.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_005.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_006.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_006.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_007.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_007.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_008.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_008.fastq.gz
GfddRAD1_005_ATCACG_L007_R1_009.fastq.gz  GfddRAD1_005_ATCACG_L007_R2_009.fastq.gz
```
Now we specify both restriction enzymes using the <code>--renz_1</code> and <code>--renz_2</code> flags along with the type combinatorial barcoding used. Here is the command: 
```
process_radtags -P -p ./raw  -b ./barcodes/barcodes -o ./samples/ \
                  -c -q -r --inline_index --renz_1 nlaIII --renz_2 mluCI
```
**The output of process_radtags**

The output of the <font color="green">process_radtags</font>  differs depending if you are processing single-end or paired-end data. In the case of single-end reads, the program will output one file per barcode into the output directory you specify. If the data do not have barcodes, then the file will retain its original name.

If you are processing _paired-end reads_, then you will get four files per barcode, two for the single-end read and two for the paired-end read. For example, given barcode ACTCG, you would see the following four files: 
```
sample_ACTCG.1.fq
sample_ACTCG.rem.1.fq
sample_ACTCG.2.fq
sample_ACTCG.rem.2.fq
```
The <font color="green">process_radtags</font>  program wants to keep the reads in phase, so that the first read in the <code>sample_XXX.1.fq</code> file is the mate of the first read in the <code>sample_XXX.2.fq</code> file. Likewise for the second pair of reads being the second record in each of the two files and so on. When one read in a pair is discarded due to low quality or a missing restriction enzyme cut site, the remaining read can’t simply be output to the <code>sample_XXX.1.fq</code> or <code>sample_XXX.2.fq</code> files as it would cause the remaining reads to fall out of phase. Instead, this read is considered a remainder read and is output into the <code>sample_XXX.rem.1.fq</code> file if the paired-end was discarded, or the <code>sample_XXX.rem.2.fq</code> file if the single-end was discarded.

**Modifying how process_radtags executes**

The <font color="green">process_radtags</font>  program can be modified in several ways. If your data do not have barcodes, omit the barcodes file and the program will not try to demultiplex the data. You can also disable the checking of the restriction enzyme cut site, or modify what types of quality are checked for. So, the program can be modified to only demultiplex and not clean, clean but not demultiplex, or some combination.

There is additional information available in [<font color="green">process_radtags</font>  manual page](https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php)

**Example: Processing data from NCBI’s Short Read Archive**

<font color="green">process_radtags</font>  can be used to process public data available in SRA and other short-read sequence databases; however, this often requires providing specific parameters when running the program. 
- SRA data is often available separately for each sample, therefore de-multiplexing may not be required. In such cases, <font color="green">process_radtags</font>  is run separately for each sample, omitting the barcodes (-b/--barcodes) option. 
- FASTQ files downloaded from SRA often do not follow Illumina’s file naming convention and should therefore be specified as generic FASTQ file inputs to process_radtags. This can be done using the -f and -1/-2 options, depending on if the downloaded data is single- or paired-end, respectively. In such cases, the --basename option can be used to provide a specific prefix name to the output files and logs. Output files and logs are otherwise named according to the specified inputs. 
- FASTQ files downloaded from SRA often do not follow Illumina’s file naming convention and should therefore be specified as generic FASTQ file inputs to process_radtags. This can be done using the -f and -1/-2 options, depending on if the downloaded data is single- or paired-end, respectively. In such cases, the --basename option can be used to provide a specific prefix name to the output files and logs. Output files and logs are otherwise named according to the specified inputs.
- By default, process_radtags will use the information available on the input file’s FASTQ headers to detect the specific sequencing platform. This information is used to identify poly-G artefacts when filtering reads. These intact FASTQ headers might not be available in FASTQ files downloaded from SRA, which will lead to poly-G runs being ignored during filtering. The user can specify the <code>--force-poly-g-check</code> option, which will result in the filtering of reads with poly-G runs, independent of the sequencing platform information. 
- Additional parameters might be needed depending on how the data was processed prior to submission to SRA. For example, the data might contain variable read lengths (which can be made uniform using <code>-t/--</code> truncate along with <code>--len-limit</code>) or the restriction enzyme cut site was previously removed, which may require disabling the cut site check (<code>--disable-rad-check</code>).

Here, we provide an example of how to process RADseq data available on SRA (accession number [SRR20082702](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR20082702&display=metadata)). The data is for an single mackarel icefish individual, paired-end sequenced on a single-digest RADseq library using the enzyme SbfI ([Rivera-Colón et al. 2023](https://doi.org/10.1093/molbev/msad029)). 

First, we download the data from NCBI using the fastq-dump program, available as part of NCBI’s [SRA Toolkit](https://github.com/ncbi/sra-tools). 

```
fastq-dump --accession SRR20082702 --gzip --split-files \
    --defline-seq '@$ac.$si/$ri' --defline-qual '+$ac.$si/$ri' --outdir ./raw
```

After downloading the paired files from this accession, named <code>SRR20082702_1.fastq.gz</code> and <code>SRR20082702_2.fastq.gz</code>, we can proceed with process_radtags. We want to apply quality filters to this data (<code>--clean</code> and <code>--quality</code>), specifying and rescueing SbfI cut site (<code>--renz-1 sbfI</code> and <code>--rescue</code>), and manually enabling the filtering of poly-G runs (<code>--force-poly-g-check</code>). Also, we want to apply a more informative name to the resulting files, which we can specify using the <code>--basename</code> option. 

```
process_radtags -1 ./raw/SRR20082702_1.fastq.gz -2 ./raw/SRR20082702_2.fastq.gz \
    --out-path ./processed --clean --quality --rescue --renz-1 sbfI \
    --force-poly-g-check --basename cgunnari_49
```


Upon completion, process_radtags will generate the usual outputs, which will be named according to the value specified by the <code>--basename</code> option. In this case, <code>cgunnari_49</code>:
```
cgunnari_49.1.fq
    cgunnari_49.rem.1.fq
    cgunnari_49.2.fq
    cgunnari_49.rem.2.fq
    cgunnari_49.log</pre>
``` 

## **Ejercicio 2: Siguiendo los pasos del tutorial anterior, genera un repositorio entro de tu cuenta de Github que se llame "Tareas_BioinfRepro2026_TusIniciales".**

Respuesta: Acorde a lo solicitado se genera el repositorio llamado [Tareas_BioinfRepro2026_JITN](https://github.com/joaquintorresnunez/Tareas_BioinfRepro2026_JITN)

## **Ejercicio 3: Clona el repositorio de la clase y actualízalo que vez que sea necesario. NOTAS IMPORTANTES PARA ESTE EJERCICIO**

```joaquin@joaquin-Victus-by-HP-Laptop-16-d0xxx:~/Escritorio$ git clone https://github.com/u-genoma/BioinfinvRepro.git --branch master --single-branch
Clonando en 'BioinfinvRepro'...
remote: Enumerating objects: 3601, done.
remote: Counting objects: 100% (260/260), done.
remote: Compressing objects: 100% (106/106), done.
remote: Total 3601 (delta 194), reused 194 (delta 151), pack-reused 3341 (from 2)
Recibiendo objetos: 100% (3601/3601), 162.45 MiB | 4.35 MiB/s, listo.
Resolviendo deltas: 100% (2034/2034), listo.
joaquin@joaquin-Victus-by-HP-Laptop-16-d0xxx:~/Escritorio$ cd BioinfinvRepro
joaquin@joaquin-Victus-by-HP-Laptop-16-d0xxx:~/Escritorio/BioinfinvRepro$ git pull
Ya está actualizado.
joaquin@joaquin-Victus-by-HP-Laptop-16-d0xxx:~/Escritorio/BioinfinvRepro$ git merge
```

## **Ejercicio Genera un repositorio dentro de tu cuenta de Github que se llame "Tareas_BioinfRepro2019_TusIniciales".**

Respuesta: Acorde a lo solicitado se genera el repositorio llamado [Tareas_BioinfRepro2019_JITN](https://github.com/joaquintorresnunez/Tareas_BioinfRepro2019_JITN)

## **Ejercicio Agrégme a mi como colaborador en el repositorio de tareas del curso que creaste en tu cuenta de Github. Mi nobre de usuario es "ravuch"**

Respuesta:
![alt text](<Captura desde 2026-08-30 16-36-25.png>)

## **Ejercicio: Mira el siguiente script (tomado del manual de Stacks) y contesta lo siguiente:**

1. ¿Cuántos pasos tiene este script?
    - Align with GSnap and convert to BAM
    - Run Stacks on the gsnap data; the i variable will be our ID for each sample we process.
    - Use a loop to create a list of files to supply to cstacks.
    - Build the catalog; the "&>>" will capture all output and append it to the Log file.
    - Calculate population statistics and export several output files.
	
    **Consta de 5 pasos**

2. ¿Si quisieras correr este script y que funcionara en tu propio equipo, qué línea deberías cambiar y a qué?
    - Reemplazar por propia dureccion a diectorio: src=$HOME/research/project 
    - Reemplazar por propio nombres de archivos files=”sample_01 sample_02 sample_03” 
    - Ruta del directorio donde está la base de datos: -A sam -d gac_gen_broads1_e64 \
            -D ~/research/gsnap/gac_gen_broads1_e64 \ 

3. ¿A qué equivale $HOME?

    Es la variable de entorno del sistema operativo (Linux) que representa la ruta absoluta al directorio personal en la terminal 

4. ¿Qué paso del análisis hace el programa gsnap?

    Realiza el alineamiento o mapeo de las lecturas de secuenciación (archivos FASTQ/FQ) contra un genoma de referencia.

5. ¿Qué hace en términos generales cada uno de los loops?

    1. Para cada archivo en el directorio files se alinea con una secuencia de referencia y lo guarda en archivo .bam/.sam
```
# Align with GSnap and convert to BAM
# 
for file in $files
do
    gsnap -t 36 -n 1 -m 5 -i 2 --min-coverage=0.90 \
            -A sam -d gac_gen_broads1_e64 \
            -D ~/research/gsnap/gac_gen_broads1_e64 \
            $src/samples/${file}.fq > $src/aligned/${file}.sam
    samtools view -b -S -o $src/aligned/${file}.bam $src/aligned/${file}.sam 
    rm $src/aligned/${file}.sam 

```

2. Inicializa una variable contador en 1, y para cada archivo en el directorio files se utiliza la función pstacks  y finalmente se mueve al archivo siguiente con i+1

```    
    # Run Stacks on the gsnap data; the i variable will be our ID for each sample we process.
# 
i=1 
for file in $files 
do 
    pstacks -p 36 -t bam -m 3 -i $i \
              -f $src/aligned/${file}.bam \
              -o $src/stacks/ 
    let "i+=1"; 
done 
```

3. Se crea una variable llamada ´samp´, y por cada archivo en files el producto previo de la función pstacks es guardado en dicha variable como una cadena de string.

``` 
 Use a loop to create a list of files to supply to cstacks.
# 
samp="" 
for file in $files 
do 
    samp+="-s $src/stacks/$file "; 
done 
```
 	
4. Itera ordenadamente por cada una de las muestras en directorio files. Utiliza la función sstacks para cada uno de los arhivos creados ‘samp’, y se genera un archivo de salidas .log

```# 
# Build the catalog; the "&>>" will capture all output and append it to the Log file.
# 
cstacks -g -p 36 -b 1 -n 1 -o $src/stacks $samp &>> $src/stacks/Log 

for file in $files 
do 
    sstacks -g -p 36 -b 1 -c $src/stacks/batch_1 \
             -s $src/stacks/${file} \ 
             -o $src/stacks/ &>> $src/stacks/Log 
done
``` 


## **Ejercicio: retoma el ejercicio anterior y divídelo en un subscript para cada paso y un script maestro que corra toda la pipeline.**

##Este script alinea las lecturas de secuenciación de cada muestra contra un genoma de referencia usando GSnap y Samtools. A continuación, procesa los alineamientos mediante Stacks para identificar loci genómicos y detectar SNPs. Calcula estadísticas de genética de poblaciones.

##Running sub-script

```
#!/bin/bash
for file in $files
do
    gsnap -t 36 -n 1 -m 5 -i 2 --min-coverage=0.90 \
            -A sam -d gac_gen_broads1_e64 \
            -D ~/research/gsnap/gac_gen_broads1_e64 \
            $src/samples/${file}.fq > $src/aligned/${file}.sam
    samtools view -b -S -o $src/aligned/${file}.bam $src/aligned/${file}.sam 
    rm $src/aligned/${file}.sam 
done
```

```
#!/bin/bash
i=1 
for file in $files 
do 
    pstacks -p 36 -t bam -m 3 -i $i \
              -f $src/aligned/${file}.bam \
              -o $src/stacks/ 
    let "i+=1" 
done
```

```
#!/bin/bash
samp="" 
for file in $files 
do 
    samp+="-s $src/stacks/$file " 
done
cstacks -g -p 36 -b 1 -n 1 -o $src/stacks $samp &>>$src/stacks/Log
```

```
#!/bin/bash
for file in $files 
do 
    sstacks -g -p 36 -b 1 -c $src/stacks/batch_1 \
             -s $src/stacks/${file} \
             -o $src/stacks/ &>>$src/stacks/Log 
done
```

```
#!/bin/bash
populations -t 36 -b 1 -P $src/stacks/ -M$src/popmap \
              -p 9 -f p_value -k -r 0.75 -s --structure --phylip --genepop
```


##Running script master

```
#!/bin/bash
export src=$HOME/research/project
export files="sample_01 sample_02 sample_03"

mkdir -p "$src/aligned" "$src/stacks"

echo "--- Iniciando Pipeline de Stacks ---"

echo "[1/5] Ejecutando GSnap)"
bash 01_align.sh

echo "[2/5] Ejecutando pstacks"
bash 02_pstacks.sh

echo "[3/5] Ejecutando cstacks"
bash 03_cstacks.sh

echo "[4/5] Ejecutando sstacks"
bash 04_sstacks.sh

echo "[5/5] Ejecutando populations"
bash 05_populations.sh

echo "Pipeline finalizada con éxito"
```