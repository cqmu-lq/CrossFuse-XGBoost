import numpy as np
import pickle

from sklearn import preprocessing,metrics
from sklearn.model_selection import GridSearchCV,train_test_split # modules for parameter tunning
from sklearn.model_selection import RepeatedKFold,cross_validate,cross_val_score # modules for cross-validation
# from lazypredict.Supervised import LazyClassifier, LazyRegressor

"""public methods"""
# Read data from a pkl file.
def read_pkl(filename):
    with open(filename, 'rb') as f:
        data=pickle.load(f)
    return data

# Write data into a pkl file.
def write_pkl(filename,data):
    try:
        with open(filename,'wb') as f:
            pickle.dump(data,f)
        return 'Succeed'
    except Exception as e:
        print('error:',e)

# ======== print statistics
# @param score_list：a list of scores
# @return: the mean value, the standard deviation, the maximum and the minimum with 4 digits after the decimal point.
def output_res(score_list,algorithm_name):
    score_array=np.array(score_list)
    avg_data=np.around(np.average(score_array),4)
    std_data=np.around(np.std(score_array,ddof=1),4)
    min_data=np.around(min(score_array),4)
    max_data=np.around(max(score_array),4)
    print(algorithm_name,"  avg:",avg_data,"std:",std_data," min:",min_data," max:",max_data)
    
# ======= Normalization methods
# @param X_data：feature vectors
# @return：standarized feature vectors
from sklearn import preprocessing
def standard_norm(X_data):
    # Standarization process
    std = preprocessing.StandardScaler()
    X_data= std.fit_transform(X_data)

    # Normalization process
    mms=preprocessing.MinMaxScaler()
#     mms=preprocessing.MaxAbsScaler()
    X_data=mms.fit_transform(X_data)
    print("The dimensionality after normalization: ",X_data.shape)
    return X_data


"""Model evaluation method"""
# ======== Cross-validation
# @param model_dict：the dictionary of models
# @param eval_score: the name of evaluation index
# @param X_all: input features 
# @param y_all: data labels
# @param n_jobs: the number of jobs to run in parallel. -1 means to use all processors.
def cross_val(model_dict,eval_score,X_all,y_all,n_jobs):
    for name,model in model_dict.items():
        cv=RepeatedKFold(n_splits=10, n_repeats=1, random_state=10)
        cross_res=cross_val_score(model, X_all, y_all, cv=cv,n_jobs=n_jobs,scoring=eval_score)
        output_res(cross_res,name)
        

# ======== Cross-validation using the correlation coefficient
# @param model_dict：the dictionary of models
# @param X_all: input features 
# @param y_all: data labels
from sklearn.model_selection import RepeatedKFold
def correlation_cross_val(model_dict,X_all,y_all):
    for name,model in model_dict.items():
        print("=======",name)
        cal_correlation(name,model,X_all,y_all)
        
from scipy import stats
def cal_correlation(reg_modelname,reg_model,X_all,y_all):
    pearsonr_corr_list=[]
    spearmanr_corr_list=[]
    train_index_list,test_index_list,Y_test_list,Y_pre_list=[],[],[],[] # stores the data set for each iteration of cross validation.
    df_2list=[]# stores the feature importance for each iteration of cross validation.
    
    cv=RepeatedKFold(n_splits=40, n_repeats=1, random_state=12)
    for (train_index, test_index) in cv.split(X_all):
        # 
        if isinstance(X_all,pd.core.frame.DataFrame):
            X_train=X_all.iloc[train_index,:]
            X_test=X_all.iloc[test_index,:]
        elif isinstance(X_all,np.ndarray):
            X_train=X_all[train_index,:]
            X_test=X_all[test_index,:]
        Y_train=y_all[train_index].ravel()
        Y_test=y_all[test_index].ravel()
        
        reg_model.fit(X_train,Y_train)
        Y_pre=reg_model.predict(X_test)
        
        pearsonr_corr=stats.pearsonr(Y_pre,Y_test)[0]
        pearsonr_corr_list.append(pearsonr_corr)
        
        spearmanr_corr=stats.spearmanr(Y_pre,Y_test)[0]
        spearmanr_corr_list.append(spearmanr_corr)
        
        # generates a 2-dimensional list with a shape of (40, n)
        train_index_list.append(train_index.tolist())
        test_index_list.append(test_index.tolist())
        Y_test_list.append(Y_test.tolist())
        Y_pre_list.append(Y_pre.tolist())
        df_2list.append(reg_model.feature_importances_)
    
    output_res(pearsonr_corr_list,'')
    print(np.round(pearsonr_corr_list,4))
    return train_index_list,test_index_list,Y_test_list,Y_pre_list,df_2list
    

"""Get standardized SMILES string"""
from rdkit import Chem
import deepchem as dc
# ======== standardized smiles
# @param X_smiles：old smiles list
# @return：new standardized smiles list
def standard_smi(X_smiles):
    X_smiles_standard=[]
    for one_smile in X_smiles:
        cs=Chem.CanonSmiles(one_smile)
        X_smiles_standard.append(cs)
    print('smiles length:',len(X_smiles_standard))
    return X_smiles_standard
# ======== Generate molecular fingerprints
# @param X_smiles_standard：standarized smiles
# @param filepath: the save directory
# @return: all kinds of molecular fingerprints
import numpy as np
def make_fps(X_smiles_standard,X_canonical_smiles,filepath):
    
    # Circular Fingerprint
    featurizer = dc.feat.CircularFingerprint(size=600)
    Xfp_circul=xsmiles_fp_list1=featurizer.featurize(X_smiles_standard)
    np.save(filepath+'CircularFingerprint.npy',Xfp_circul)
    print("CircularFingerprint：",Xfp_circul.shape)

    # MACCS Keys Fingerprint
    featurizer = dc.feat.MACCSKeysFingerprint()
    Xfp_maccsk=featurizer.featurize(X_smiles_standard)
    np.save(filepath+'MACCSKeysFingerprint.npy',Xfp_maccsk)
    print("MACCSKeysFingerprint：",Xfp_maccsk.shape)

    # PubChem Fingerprint--881 columns
    featurizer = dc.feat.PubChemFingerprint()
    Xfp_pubchem = featurizer.featurize(X_canonical_smiles)
    np.save(filepath+'PubChemFingerprint.npy',Xfp_pubchem)
    print("PubChemFingerprint：",Xfp_pubchem.shape)

    #RDKit Descriptors
    featurizer = dc.feat.RDKitDescriptors()
    Xfp_rdkit=featurizer.featurize(X_smiles_standard)
    np.save(filepath+'RDKitDescriptors.npy',Xfp_rdkit)
    print("RDKitDescriptors：",Xfp_rdkit.shape)
    
    # Mordred Descriptors
    featurizer = dc.feat.MordredDescriptors()
    Xfp_md_3d = featurizer.featurize(X_smiles_standard)
    featurizer = dc.feat.MordredDescriptors(ignore_3D=True)
    Xfp_md = featurizer.featurize(X_smiles_standard)
    np.save(filepath+'MordredDescriptors.npy',Xfp_md_3d)
    print("MordredDescriptors：",Xfp_md_3d.shape)
    
    return Xfp_circul,Xfp_maccsk,Xfp_pubchem,Xfp_rdkit,Xfp_md_3d,Xfp_md

# ======== Transform molecular fingerprints to dataframe. Adjust its column name and save it in csv format.
# @param fp_name：the column name
# @param fp_array: the array of molecular fingerprints
# @param filepath: the save directory
# @return：a dataframe of the given molecular fingerprints
import pandas as pd
def array2df(fp_name,fp_array,filepath):
    row_num,col_num=fp_array.shape
    col_name_list=[]
    for i in range(col_num):
        col_name_list.append(fp_name+'_'+str(i))
    fp_df=pd.DataFrame(fp_array)
    fp_df.columns=col_name_list
    filename=fp_name+'.csv'
    fp_df.to_csv(filepath+filename)
    print(filepath+filename)
    return fp_df

"""regression model"""
from sklearn.linear_model import LinearRegression,Lars # linear regression
from sklearn.neighbors import KNeighborsRegressor # k-nearest neighbor regression
from sklearn.tree import DecisionTreeRegressor # regression tree

from sklearn.svm import SVR # support vector machine regression
from sklearn.neural_network import MLPRegressor # ANN
# Ensemble learning: Adaboost, Random Forest, Extremely Random Forest, Gradient Boosting Tree
from sklearn.ensemble import AdaBoostRegressor,RandomForestRegressor,ExtraTreesRegressor,GradientBoostingRegressor,BaggingRegressor
from xgboost import XGBRegressor #XGBoost

def regression_model(n_jobs):
    model_dict={}
# Linear regression
#     lr_model=LinearRegression(n_jobs=6) 
#     model_dict["LinearRegressor"]=lr_model
# Ridge regression
#     lars_model=Lars() 
#     model_dict["RidgeRegressor"]=lars_model
    # k-nearest neighbor regression
    knr_model=KNeighborsRegressor(weights="uniform",n_jobs=13) 
    model_dict["KNN"]=knr_model
    # Regression Tree
    dtreer_model=DecisionTreeRegressor() 
    model_dict["RegTree"]=dtreer_model
# SVM regression
#     svr_model=SVR(kernel="linear") 
#     model_dict["SVR"]=svr_model
    # ANN
    mlpr_model=MLPRegressor(max_iter=5000) 
    model_dict["ANN"]=mlpr_model
    
# Adaboost
    adaboost_model=AdaBoostRegressor() 
    model_dict["AdaBoost"]=adaboost_model
# Random forest
    randomfr_model=RandomForestRegressor(n_jobs=n_jobs,random_state=10) 
    model_dict["RF"]=randomfr_model
# Extremely Random Forest
    etreer_model=ExtraTreesRegressor(n_jobs=n_jobs,random_state=10) 
    model_dict["ExtraTrees"]=etreer_model
# Gradient boosting tree
    gbr_model=GradientBoostingRegressor() 
    model_dict["GBoost"]=gbr_model
# XGBoost
    xgboostr_model=XGBRegressor(n_jobs=n_jobs,random_state=10) 
    model_dict["XGBoost"]=xgboostr_model
# Bagging
#     bgr_model=BaggingRegressor(n_jobs=n_jobs,random_state=10) 
#     model_dict["BaggingRegressor"]=bgr_model

    return model_dict

"""classification model"""
from sklearn.linear_model import LogisticRegression # Logistic
from sklearn.neighbors import KNeighborsClassifier # KNN
from sklearn.naive_bayes import MultinomialNB # Naive Bayes
from sklearn.tree import DecisionTreeClassifier # Decision tree
from sklearn.svm import SVC # SVM
from sklearn.neural_network import MLPClassifier # ANN
# Ensemble learning: Adaboost, Random Forest, Extremely Random Forest, Gradient Boosting Tree
from sklearn.ensemble import AdaBoostClassifier,RandomForestClassifier,ExtraTreesClassifier,GradientBoostingClassifier,BaggingClassifier
from xgboost import XGBClassifier #XGBoost
def classification_model(n_jobs):
    model_dict={}

#     lr_model=LogisticRegression(max_iter=1000) # Logistic
#     model_dict["Logistic"]=lr_model
    
#     knn_model=KNeighborsClassifier() # KNN
#     model_dict["KNN"]=knn_model
    
#     bys_model=MultinomialNB() # Naive Bayes
#     model_dict["NB"]=bys_model
    
#     dtree_model=DecisionTreeClassifier() # Decision tree
#     model_dict["DecisionTree"]=dtree_model
    
#     svm_model=SVC() # SVM
#     model_dict["SVM"]=svm_model
    
#     mlpc_model=MLPClassifier(max_iter=3000) # ANN
#     model_dict["ANN"]=mlpc_model

    adaboostc_model=AdaBoostClassifier(random_state=10)
    model_dict["AdaBoost"]=adaboostc_model

    randomfc_model=RandomForestClassifier(n_jobs=n_jobs,random_state=10)
    model_dict["RF"]=randomfc_model
    
    etreec_model=ExtraTreesClassifier(random_state=10) # Extremely Random Forest
    model_dict["ExtraTrees"]=etreec_model
    
    gbc_model=GradientBoostingClassifier(random_state=10) # Gradient boosting tree
    model_dict["GBoost"]=gbc_model
    
    xgboostc_model=XGBClassifier(n_jobs=n_jobs,random_state=10) #XGBoost
    model_dict["XGBoost"]=xgboostc_model
    
    bgc_model=BaggingClassifier(n_jobs=n_jobs,random_state=10) #Bagging
    model_dict["Bagging"]=bgc_model
    
    return model_dict
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
