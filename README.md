# Addressing Informative Presence in Electronic Health Records Data when Implementing Prognostic Models

This repository contains the source code, simulation frameworks, and supplementary materials for the paper **"Addressing Informative Presence in Electronic Health Records Data when Implementing Prognostic Models"**. 

The project evaluates the statistical implications of informative presence—specifically focusing on the bias introduced by traditional imputation strategies (such as fixed-value imputation)—and demonstrates robust methods for prediction using Electronic Health Records (EHR) data. The application focus of the main simulation and validation workflow centers around prognostic modeling for hospital-acquired Venous Thromboembolism (HA-VTE) risk.

---

## Repository Structure & File Descriptions

### 1. `Proof.pdf`
* **Description:** Hand-written analytical proof providing the theoretical derivation of statistical bias when using fixed-value imputation under different missing mechanisms.
* **Role:** Provides additional technical and mathematical details behind the bias derivations, laying the foundational theory the mathematical foundation for the expected direction and magnitude of estimation bias, justifying the empirical results observed in the subsequent simulation studies.

### 2. `bias_simulation_missing.R`
* **Description:** R script dedicated to executing the simulation study that empirically validates the analytical bias derivation presented in `Proof.pdf`.
* **Role:** This script systematically generates synthetic data under varying degrees of missingness and informative presence mechanisms, applies fixed-value imputation, and estimates the resulting bias to confirm that the theoretical derivations align with empirical data.

### 3. `function.R`
* **Description:** A centralized R script containing all custom helper functions, data-generating processes, metrics calculations, and model execution routines utilized throughout the simulation framework.
* **Role:** This file acts as the primary functional backend. Sourcing this file ensures reproducible execution and code cleanliness across both the bias simulation script and the main simulation workflows.

### 4. `vte_simulation.Rmd`
* **Description:** The primary reproducible R Markdown notebook covering the end-to-end main simulation study for the paper.
* **Key Components:**
  * **Data Generation:** Randomly samples covariate information and simulates outcomes and absences based on outcome-generating models and absence-generating models.
  * **Imputation Strategy:** Implements and compares multiple conventional imputation approaches.
  * **Prediction Modeling:** Trains prognostic models to predict Venous Thromboembolism (VTE) risk.
  * **Result Visualization:** Generates the core tables and figures included in the results summary section of the manuscript.

### 5. `vte_validation.Rmd`
* **Description:** R Markdown notebook containing the data processing and statistical workflows for the temporal validation cohort analysis.
* **Role:** Corresponds directly to the **Application** section of the paper. It tests the transportability, calibration, and discrimination of the developed VTE prognostic models in an independent validation data structure derived from clinical EHR environments, accounting for the informative presence of laboratory results.

---

## Getting Started

### Prerequisites
To run these analyses, you will need an active installation of **R** along with standard packages for data manipulation, predictive modeling, and Markdown reporting. 

Required packages typically include:
```R
install.packages(c("tidyverse", "rmarkdown", "knitr", "survival", "pROC", "ResourceSelection"))
# Additional machine learning or evaluation packages used within function.R
