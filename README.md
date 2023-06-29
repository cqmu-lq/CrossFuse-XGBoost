# Data Availability of CrossFuse-XGBoost  
### CrossFuse-XGBoost: Accurate prediction of maximum daily recommended dose through multi-feature fusion and cross-validation screening  
### Preamble

This repository contains the data and code for the paper titled "CrossFuse-XGBoost: Accurate prediction of maximum daily recommended dose through multi-feature fusion and cross-validation screening".

## Abstract
In the process of drug development, approximately 30% of failures are attributed to drug safety issues. Among them, the first-in-human (FIH) trial of a new drug represents one of the highest safety risks, and the selection of the initial dose is crucial for ensuring the safety of clinical trials. Traditional dose estimation methods, when extrapolating from animal data to humans, have resulted in catastrophic events during Phase I clinical trials due to interspecies differences in compound sensitivity and unknown molecular mechanisms.

To address this issue, this study proposes a CrossFuse-XGBoost method that can directly predict the maximum daily recommended dose of a compound based on existing human research data, providing a reference for FIH dose selection. This method not only integrates multiple features, including molecular representations, physicochemical properties, and compound-protein interactions but also improves feature selection based on cross-validation.

The results demonstrate that compared to existing local weighted methods (k-nn and v-nn), the CrossFuse-XGBoost method not only improves prediction accuracy but also solves the low prediction coverage issue of v-nn, achieving full coverage of the test set and enabling more reliable predictions. Furthermore, this study offers a high level of interpretability by identifying the importance of different features in model construction. A total of 241 features that have the most significant impact on the maximum daily recommended dose were selected, providing references for optimizing the structure of new compounds and guiding experimental research.

## Data Availability
The data used in this study, along with the software code, are made available in this repository to ensure transparency and reproducibility of the findings. The following files are included:

1. **[dataset.csv](link_to_dataset.csv)**: The dataset used for training and evaluation. Each row represents a compound with its associated features and the maximum daily recommended dose.

2. **[crossfuse_xgboost.py](link_to_code.py)**: The Python code for implementing the CrossFuse-XGBoost method and performing the prediction.

## Usage
To use the CrossFuse-XGBoost method and reproduce the results, follow these steps:

1. Clone or download this repository to your local machine.

2. Install the required dependencies listed in [requirements.txt](link_to_requirements.txt) using the package manager of your choice.

3. Run the [crossfuse_xgboost.py](link_to_code.py) script, providing the path to the dataset as an input. Make sure the dataset is in the appropriate format as specified in the script.

4. The script will execute the CrossFuse-XGBoost method and generate the predicted maximum daily recommended dose for each compound in the dataset.

Please refer to the code comments for further instructions on customizing the method or adapting it to your specific use case.

## Citation
If you find this work useful or build upon it, please consider citing:


## License
[Specify the license you choose for your code and data]

## Contact
For any inquiries or questions regarding this research or the code, please contact Jianbo Pan at panjianbo@cqmu.edu.cn.

---
Note: Replace the placeholders [link_to_dataset.csv], [link_to_code.py], [link_to_requirements.txt], [Your Name], and [Your Email Address] with the actual links and contact information. Also, ensure you choose an appropriate license for your code and data,

