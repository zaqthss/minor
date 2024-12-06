# MINOR

Source code and supplementary materials for "MINOR: Multivariate Time Series Iterative Cleaning Algorithm".

To facilitate reproducibility, we have open-sourced all the datasets, algorithms, and code presented in the paper, and this document provides a guideline for reproduction.

We used publicly available code sources for the deep leaning baselines: 

* USAD: https://github.com/imperial-qore/TranAD
* TranAD: https://github.com/imperial-qore/TranAD

## **Note**

* Code for all the algorithms in the paper (except for the deep learning baselines) are located in the `MINOR` directory. 
* The code in the `MINOR_py` directory is used for the application experiments in subsection D of Section V.
* Please open these two directories using Java and Python IDEs, respectively.

## Requirement

* java
  * java : 17
  * IDE : IntelliJ IDEA
  * Maven is required to download the libraries required by the code.
* python (IDE : pycharm)
  * python==3.10
  * pandas==2.2.2
  * scikit-learn==1.5.1

## Code Structure

* MINOR/
  * data/: This directory stores all the datasets used in this paper. 
    * app/: This directory stores datasets with errors for application experiments in Section V.D
    * dirty/: This directory stores datasets with errors
    * raw/: 
      * Classification/: This directory stores original datasets for classification experiments in Section V.D
      * Cluster/: This directory stores original datasets for cluster experiments in Section V.D
  * result/: This directory stores all the result of experiments in this paper, and it can be reproduced according to the following "Experiments". 
  * src/com/MTCSC/:
    * core/:
      * AkaneHeuristic.java: Corresponding to the baseline Akane
      * AR.java: Implementation of the AR model
      * ARX.java: Implementation of the ARX model
      * BaseCleaningModel.java: Parent class, provide some basic properties and methods
      * IMR.java: Implementation of the IMR model
      * IMRStream.java: Implementation of the IMR-stream model for online computing
      * MINOR_B.java: Implementation of the MINOR-B model
      * MINOR_base.java: Child class of MTSCModel, defined basic properties and methods for MINOR-B, MINOR-U and MINOR-O
      * MINOR_BUni.java: Implementation of the MINOR-BUni model
      * MINOR_O.java: Implementation of the MINOR-O model for online computing
      * MINOR_U.java: Implementation of the MINOR-U model
      * MTCSC.java:  Corresponding to the baseline MTCSC-C
      * MTCSC_A.java:  Corresponding to the baseline MTCSC-A
      * MTCSC_Uni.java:  Corresponding to the baseline MTCSC-Uni for application experiments
      * MTSCModel.java: Child class of BaseCleaningModel, defined basic properties and methods for MINOR-B, MINOR-U, MINOR-O, VAR, VARX and VARX-stream
      * UTSCModel.java: Child class of BaseCleaningModel, defined basic properties and methods for AR, ARX, IMR, IMR-stream and MINOR-BUni
      * VAR.java: Implementation of the VAR model
      * VARX.java: Implementation of the VARX model
      * VARXStream.java: Implementation of the VARX-stream model for online computing
    * entity/: Defined the required classes
    * enums/: Defined the required enums
    * example/: Corresponding to the examples from example 2 to example 9
    * experiment/: The directory stores all the Java files corresponding to the experiments in this paper. Details can be found in the following "Experiments".
    * pre/: Classes of generating synthetic errors and generating random labels
    * Utils/ : Defined the required utils.

## Experiments

**Note: Before reproduce results, src/com/MTCSC/pre/RunAllInjection.java is required to be runed first to generate datasets with errors and labels.**

To reproduce the experimental results in this paper, just run a separate Java file. The Java file corresponds to each experimental section as follows.

**Path = MINOR/src/com/MTCSC/experiment/**

**Data = MINOR/data/**

**Result = MINOR/result/**

* Section V.B.1 Case Study
  * Path/`CaseStudy_Figure.java`, results in Result/CaseStudy/
  * Path/`CaseStudy_Table.java`, results in Result/CaseStudy/table/
* Section V.B.2 Varying Order p
  * Path/`P_GPS.java`, results in Result/GPS/summary/p/
* Section V.B.3 Varying Convergence Threshold τ
  * Path/`Threshold_GPS.java`, results in Result/GPS/summary/threshold/
* Section V.B.4 Specifying Maximum Number of Iterations
  * Path/`MaxIteration_GPS.java`, results in Result/GPS/summary/maxiteration/
* Section V.B.5 Varying Labeling Rate
  * Path/`LabelRate_GPS.java`, results in Result/GPS/summary/lr/
* Section V.C.1 Varying Error Rate e%
  * Path/`ErrorRate_ILD.java`, results in Result/ILD/summary/er/
  * Path/`ErrorRate_ECG.java`, results in Result/ECG/summary/er/
* Section V.C.2 Varying Error Length
  * Path/`ErrorLength_ILD.java`, results in Result/ILD/summary/el/
  * Path/`ErrorLength_ECG.java`, results in Result/ECG/summary/el/
* Section V.C.3 Varying Data Size n
  * Path/`DataSize_ILD.java`, results in Result/ILD/summary/size/
  * Path/`DataSize_ECG.java`, results in Result/ECG/summary/size/
* Section V.C.4 Varying Error Pattern
  * Path/`ErrorPattern_ILD.java`, results in Result/ILD/summary/pattern/
* Section V.C.5  Evaluation on Online Computing
  * Path/`OnlineComputing.java`, results in Result/online/
* Section V.D Applications: The final classification and clustering experiments were completed using Python. First, you need to run the corresponding Java files below to generate repair results, and then copy the folder named `repaired` in `Data/app/` to MINOR_py/data/ (already created)
  * java
    * Path/Application.java, repair results in Data/app/repaired/
  * python
    * MINOR_py/Application/classification.py:/ reproduce the results of classification experiments on `Car` and `Lightning2`
    * MINOR_py/Application/cluster.py:/ reproduce the results of cluster experiments on `DistalPhalanxTW` and `InsectEPGRegularTrain`