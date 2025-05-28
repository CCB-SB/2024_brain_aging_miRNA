### Code for the figures of the publication:

# A spatio-temporal brain miRNA expression atlas identifies sex-independent age-related microglial driven miR-155-5p increase

Annika Engel<sup>&,1</sup>, Viktoria Wagner<sup>&,1,2</sup>, Oliver Hahn<sup>2,3</sup>, Aulden G. Foltz<sup>2</sup>, Micaiah Atkins<sup>2</sup>, Amila Beganovic<sup>1</sup>, Ian H. Guldner<sup>2,4</sup>, Nannan Lu<sup>2,4</sup>, Aryaman Saksena<sup>2</sup>, Ulrike Fischer<sup>5</sup>, Nicole Ludwig<sup>5</sup>, Eckart Meese<sup>5</sup>, Tony Wyss-Coray<sup>2,4,6</sup>, Andreas Keller<sup>1,7</sup>

<sup>1</sup> Clinical Bioinformatics, Saarland University, 66123, Saarbrücken, Germany.

<sup>2</sup> Department of Neurology and Neurological Sciences, Stanford University, Stanford, CA, 94305, USA.

<sup>3</sup> Calico Life Sciences LLC, San Francisco, CA, USA.

<sup>4</sup> Wu Tsai Neurosciences Institute, Stanford University School of Medicine, Stanford, CA, USA

<sup>5</sup> Department of Human Genetics, Saarland University, 66421 Homburg/Saar, Germany.

<sup>6</sup> The Phil and Penny Knight Initiative for Brain Resilience, Stanford University, Stanford, CA, USA.

<sup>7</sup> Helmholtz Institute for Pharmaceutical Research Saarland, Helmholtz Center for Infection Research, 66123, Saarbrücken, Germany.

<sup>&</sup> These authors contributed equally: Annika Engel, Viktoria Wagner

PMID: 40382330 PMCID: [PMC12085673](https://pmc.ncbi.nlm.nih.gov/articles/PMC12085673/) DOI: [10.1038/s41467-025-59860-6](https://doi.org/10.1038/s41467-025-59860-6)

## Usage
1. Download raw count files from GEO and put them into the "raw_data" folder as described in the README files there.

2.
    1. Run the bash script `HUMAN=0 MRNA=0 bash run.sh`.
    2. If you have access to the human data from Jager et al., Nature Scientifc data (https://doi.org/10.1038/sdata.2018.142) and Bannett et al., Sage Journals (https://doi.org/10.3233/JAD-179939). You can start the bash script "HUMAN=1 MRNA=0 bash run.sh" to also get the figures of the human data.
    3. If you have included the mRNA data from Hahn et al., Cell (https://doi.org/10.1016/j.cell.2023.07.027). You can start the bash script "HUMAN=0 MRNA=1 bash run.sh" to also get the figures with the mRNA data.
    4. It is possible to combine 2.2 and 2.3 to start for all figures with `HUMAN=1 MRNA=1 bash run.sh`.

3. In the `publication` folder you find the single figures of our main and supplementary figures.
