# ECG_PLATFORM
ECG_PLATFORM is a complete framework designed for testing QRS detectors and ECG delineators on publicly available datasets.


In ECG_PLATFORM there are included the following algorithms for the detection of QRS complexes:

1) engzee segmenter - Andr ́e Louren ̧co, Hugo Silva, Paulo Leite, Renato Louren ̧co, and Ana Fred.Real  Time  Electrocardiogram  Segmentation  for  Finger  Based  ECG  Bio-metrics.
2) hamilton segmenter - Hamilton, Tompkins, W. J., "Quantitative investigation of QRS detection rules using the MIT/BIH  arrhythmia  database",  IEEE  Trans.  Biomed.  Eng.,  BME-33,  pp.  1158-1165,  1987.
3) WT delineator - Juan Pablo Martínez, Rute Almeida, Salvador Olmos, Member, IEEE, Ana Paula Rocha, and Pablo
Laguna, Member, IEEE A Wavelet-Based ECG Delineator: Evaluation on Standard Databases
4) ECGPUWAVE - https://www.physionet.org/content/ecgpuwave/1.3.4/

Data sets on which the methods are tested:

1) cinc1 - https://physionet.org/pn3/challenge/2014/set-p/
2) cinc2 - https://physionet.org/pn3/challenge/2014/set-p2/
3) mitdb - https://archive.physionet.org/physiobank/database/mitdb/
4) mitdb (pwave annotations) - https://archive.physionet.org/physiobank/database/pwave/
5) qtdb - https://physionet.org/content/qtdb/1.0.0/
6) ludb - https://physionet.org/content/ludb/1.0.0/
7) telehealth - https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/QTG0EP

Required Python packages: biosppy, wfdb, numpy, pandas, h5py, scipy, matplotlib, scikit-learn
