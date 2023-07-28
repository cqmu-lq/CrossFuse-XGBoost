# Data Availability of CrossFuse-XGBoost   
### CrossFuse-XGBoost: Accurate prediction of the maximum daily recommended dose through multifeature fusion, cross-validation screening and extreme gradient boosting
### Preamble

This repository contains the data and code for the paper titled "CrossFuse-XGBoost: Accurate prediction of the maximum daily recommended dose through multifeature fusion, cross-validation screening and extreme gradient boosting".

<img src="https://github.com/cqmu-lq/CrossFuse-XGBoost/blob/main/src/img/Figure%201.jpg" alt="CrossFuse"/><br/>


## Data Availability
The data used in this study, along with the software code, are made available in this repository to ensure transparency and reproducibility of the findings. The following files are included:

1. **src/data**: The dataset used for training and external validation sets. Each row represents a compound with its associated features and the maximum daily recommended dose.

2. **src/code**: The Python code for implementing the CrossFuse-XGBoost method and performing the prediction.

## Usage
To use the CrossFuse-XGBoost method and reproduce the results, follow these steps:

1. Clone or download this repository to your local machine.

2. Install the required dependencies listed in "src/requirements.txt" using the package manager of your choice.

3. Run every cell of "src/data/crossfuse_xgboost.ipynb" according to the instructions in Jupyter. Make sure that the data set is stored in the same path as specified in the script.

4. The script will execute the CrossFuse-XGBoost method and generate the predicted maximum daily recommended dose for each compound in the dataset. Each cell has a functional comment about the result.

Please refer to the code comments for further instructions on customizing the method or adapting it to your specific use case.

## Citation
If you find this work useful or build upon it, please consider cite it.


## Contact
For any inquiries or questions regarding this research or the code, please contact Jianbo Pan at panjianbo@cqmu.edu.cn.


## Copyright License
Permission is hereby granted, free of charge, to any person obtaining a copy of this project and associated documentation files (the "project"), to deal in the project without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the project, and to permit persons to whom the project is furnished to do so, subject to the following conditions:

The data of this study, including the dataset and associated findings, are intended for experimental reference only and should not be directly interpreted or used to guide human medicine. The authors do not assume any responsibility or liability for any consequences arising from the use or application of the data. The project is provided "as is," without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, and non-infringement. in no event shall the authors or copyright holders be liable for any claim, damages, or other liability, whether in an action of contract, tort, or otherwise, arising from, out of, or in connection with the project or the use or other dealings in the project.

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the project.

