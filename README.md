# Charting pituitary structures in neurological diseases with normative modeling

## Table of contents

1. [Introduction](#Introduction)  
2. [Datasets](#Datasets)  
3. [Models](#Models)  
   3.1 [nnU-Net automated pituitary segmentation](#Segmentation)  
   3.2 [GAMLSS normative modeling](#Normative)  
4. [Scripts](#Scripts)  
   4.1 [Required R packages and installation](#Check-and-install-packagesR)  
   4.2 [Normative curve estimation and peak age determination](#Normative-model-fitR)  
   4.3 [Bootstrap analysis of normative curves](#Bootstrap-normative-mdoel-fitR)  
   4.4 [ICV adjustment and model comparison](#ICVR)  
   4.5 [Model calibration using new datasets](#Calibration-normative-model-using-new-datasetR)  
   4.6 [Applying normative models to disease datasets](#Disease-application-normative-modelR)  
   4.7 [Applying normative models to individual-level data](#Individual-application-normative-modelR)  
   4.8 [Statistical analyses of deviation scores across diseases](#Stastical-analysis-deviations-across-diseasesR)  
   4.9 [Clinical downstream tasks](#Clinical-tasks)  
5. [License](#License)

---

## 1. Introduction

This repository provides a complete computational workflow for constructing **pituitary subregional normative models** and applying them to clinical cohorts spanning multiple neurological diseases.

The workflow consists of two major components:

1. **nnU-Net–based automated pituitary segmentation** from 3D T1-weighted MRI to extract volumes for whole pituitary and subregions (anterior pituitary, posterior pituitary, pituitary stalk).  
2. **GAMLSS-based lifespan normative modeling** for estimating centile trajectories and deriving deviation scores used for disease association analyses and prognostic modeling.

The repository enables reproducible normative modeling pipelines supporting:

- disease classification vs. matched controls  
- cognitive and clinical association analyses  
- longitudinal prognosis and survival modeling  

Example datasets, scripts, and model interfaces are provided for demonstration and benchmarking.  

---

## 2. Datasets

Example datasets are provided solely for demonstration and code testing.  
They include synthetic MRI-derived pituitary volumes and demographic variables in the same format required by the modeling workflow.  

Full multi-center datasets underlying published work cannot be publicly distributed due to institutional and privacy constraints. Requests for collaboration or normative model access should be directed to the corresponding author.

---

## 3. Models

### 3.1 nnU-Net automated pituitary segmentation  <a name="Segmentation"></a>

This repository includes code for nnU-Net–based automated segmentation of pituitary subregions from 3D T1-weighted MRI.

The segmentation pipeline:

- adopts 3D nnU-Net architecture  
- uses manually annotated pituitary masks as training labels  
- outputs subject-specific segmentation masks and volumetric measurements for:  
  - whole pituitary  
  - anterior pituitary  
  - posterior pituitary  
  - pituitary stalk  

Notes:

- model weights and training data are not released due to data-use restrictions  
- inference scripts and expected I/O formats are provided for reproducibility within controlled datasets  

### 3.2 GAMLSS normative modeling  <a name="Normative"></a>

Pre-trained normative models for pituitary subregional volumes are stored as `.rds` files.  
These models adjust centile estimation for:

- age  
- sex  
- site/scanner effects  

The normative model outputs include centile scores and deviation metrics suitable for clinical and analytical applications.

Model files are not distributed publicly; access may be requested from the authors.

---

## 4. Scripts

This repository contains a complete, modular workflow for pituitary normative modeling and clinical applications.

### 4.1 Required R packages and installation  

`Check-and-install-packages.R`

Automatically installs all required R packages for normative modeling and statistical analyses.

---

### 4.2 Normative curve estimation and peak age determination  

`Normative-model-fit.R`

Fits GAMLSS models for whole gland and subregions and derives normative lifespan trajectories and peak ages.

---

### 4.3 Bootstrap analysis of normative curves  

`Bootstrap-normative-mdoel-fit.R`

Performs bootstrap resampling to obtain uncertainty bounds and confidence intervals for estimated normative curves and peak parameters.

---

### 4.4 ICV adjustment and model comparison  

`ICV.R`

Implements ICV-adjustment methods and compares adjusted and unadjusted normative model performance and trajectory consistency.

---

### 4.5 Model calibration using new datasets  

`Calibration-normative-model-using-new-dataset.R`

Performs recalibration procedures to adapt normative models to new dataset distributions, sites, or scanner profiles.

---

### 4.6 Applying normative models to disease datasets  

`Disease-application-normative-model.R`

Applies pre-trained normative models to disease cohort MRI measurements and computes deviation scores.

---

### 4.7 Applying normative models to individual-level data  

`Individual-application-normative-model.R`

Computes deviation scores for individual patients, enabling personalized normative evaluation.

---

### 4.8 Statistical analyses of deviation scores across diseases  

`Stastical-analysis-deviations-across-diseases.R`

Performs group comparisons, effect size estimation, permutation testing, and multiple-comparison correction.

---

### 4.9 Clinical downstream tasks  

This section may include scripts for:  

- ROC-based classification using deviation features  
- cognitive association analyses  
- survival/prognostic modeling based on deviation scores  

---

## 5. License

**MIT License**

Copyright (c) 2025

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. file for full terms.
