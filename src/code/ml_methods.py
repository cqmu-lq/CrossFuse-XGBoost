import numpy as np
import pickle

"""数据处理"""
from sklearn import preprocessing,metrics #数据处理
from sklearn.model_selection import GridSearchCV,train_test_split #调参
from sklearn.model_selection import RepeatedKFold,cross_validate,cross_val_score #交叉验证
# from lazypredict.Supervised import LazyClassifier, LazyRegressor

"""公共方法"""
# 读取pkl文件
def read_pkl(filename):
    with open(filename, 'rb') as f:
        data=pickle.load(f)
    return data

# 写入pkl文件
def write_pkl(filename,data):
    try:
        with open(filename,'wb') as f:
            pickle.dump(data,f)
        return 'Succeed'
    except Exception as e:
        print('error:',e)

# ========打印输出
# 输入：传入list，输入均值，最大值，最小值，保留三位有效数字
def output_res(score_list,algorithm_name):
    score_array=np.array(score_list)
    avg_data=np.around(np.average(score_array),4)
    std_data=np.around(np.std(score_array,ddof=1),4)
    min_data=np.around(min(score_array),4)
    max_data=np.around(max(score_array),4)
    print(algorithm_name,"  avg:",avg_data,"std:",std_data," min:",min_data," max:",max_data)
    
# ========标准化和归一化方法
# 输入：需要处理的数据
# 输出：处理完的数据
from sklearn import preprocessing
def standard_norm(X_data):
    # 标准化
    std = preprocessing.StandardScaler()
    X_data= std.fit_transform(X_data)

    # 归一化处理
    mms=preprocessing.MinMaxScaler()
#     mms=preprocessing.MaxAbsScaler()
    X_data=mms.fit_transform(X_data)
    print("标准化后维度：",X_data.shape)
    return X_data


"""模型评价方法"""
# ========交叉验证
# 输入：模型字典，评价分数类型，X数据，Y数据,n_jobs
# 输出：
def cross_val(model_dict,eval_score,X_all,y_all,n_jobs):
    for name,model in model_dict.items():
        cv=RepeatedKFold(n_splits=10, n_repeats=1, random_state=10)
        cross_res=cross_val_score(model, X_all, y_all, cv=cv,n_jobs=n_jobs,scoring=eval_score)
        output_res(cross_res,name)
        

# ========交叉验证---计算相关系数
# 输入：模型字典，X数据，Y数据
# 输出：
from sklearn.model_selection import RepeatedKFold
def correlation_cross_val(model_dict,X_all,y_all):
    for name,model in model_dict.items():
        print("=======",name)
        cal_correlation(name,model,X_all,y_all)
        
from scipy import stats
def cal_correlation(reg_modelname,reg_model,X_all,y_all):
    pearsonr_corr_list=[]
    spearmanr_corr_list=[]
    train_index_list,test_index_list,Y_test_list,Y_pre_list=[],[],[],[] #存储每折交叉验证的值
    df_2list=[]# 存放每折交叉验证的特征重要性df
    
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
        
        # 生成的是（10,n）的二维列表[[],[]]
        train_index_list.append(train_index.tolist())
        test_index_list.append(test_index.tolist())
        Y_test_list.append(Y_test.tolist())
        Y_pre_list.append(Y_pre.tolist())
        df_2list.append(reg_model.feature_importances_)
    
    output_res(pearsonr_corr_list,'')
    print(np.round(pearsonr_corr_list,4))
    return train_index_list,test_index_list,Y_test_list,Y_pre_list,df_2list
    

"""获取指纹"""
from rdkit import Chem
import deepchem as dc
# ========标准化smiles
# 输入：原smiles列表
# 输出：新smiles列表
def standard_smi(X_smiles):
    X_smiles_standard=[]
    for one_smile in X_smiles:
        cs=Chem.CanonSmiles(one_smile)
        X_smiles_standard.append(cs)
    print('smiles长度:',len(X_smiles_standard))
    return X_smiles_standard
# ========生成指纹
# 输入：标准化后的smiles,要保存的路径
# 输出：不同的分子指纹
import numpy as np
def make_fps(X_smiles_standard,X_canonical_smiles,filepath):
    
    # 循环指纹
    featurizer = dc.feat.CircularFingerprint(size=600)
    Xfp_circul=xsmiles_fp_list1=featurizer.featurize(X_smiles_standard)
    np.save(filepath+'循环指纹.npy',Xfp_circul)
    print("循环指纹：",Xfp_circul.shape)

    # MACCSKeysFingerprint
    featurizer = dc.feat.MACCSKeysFingerprint()
    Xfp_maccsk=featurizer.featurize(X_smiles_standard)
    np.save(filepath+'MACCSKeysFingerprint.npy',Xfp_maccsk)
    print("MACCSKeysFingerprint：",Xfp_maccsk.shape)

    # PubChemFingerprint--881列
    featurizer = dc.feat.PubChemFingerprint()
    Xfp_pubchem = featurizer.featurize(X_canonical_smiles)
    np.save(filepath+'PubChemFingerprint.npy',Xfp_pubchem)
    print("PubChemFingerprint：",Xfp_pubchem.shape)

    #RDKitDescriptors
    featurizer = dc.feat.RDKitDescriptors()
    Xfp_rdkit=featurizer.featurize(X_smiles_standard)
    np.save(filepath+'RDKitDescriptors.npy',Xfp_rdkit)
    print("RDKitDescriptors：",Xfp_rdkit.shape)
    
    # MordredDescriptors
    featurizer = dc.feat.MordredDescriptors()
    Xfp_md_3d = featurizer.featurize(X_smiles_standard)
    featurizer = dc.feat.MordredDescriptors(ignore_3D=True)
    Xfp_md = featurizer.featurize(X_smiles_standard)
    np.save(filepath+'MordredDescriptors.npy',Xfp_md_3d)
    print("MordredDescriptors：",Xfp_md_3d.shape)
    
    return Xfp_circul,Xfp_maccsk,Xfp_pubchem,Xfp_rdkit,Xfp_md_3d,Xfp_md

# ========指纹数据转换为dataframe，修改列名，并保存为csv
# 输入：列名前缀，指纹数据的array，保存路径
# 输出：指纹数据的dataframe
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

"""回归模型"""
from sklearn.linear_model import LinearRegression,Lars #线性回归
from sklearn.neighbors import KNeighborsRegressor #k近邻回归
from sklearn.tree import DecisionTreeRegressor #回归树
from sklearn.svm import SVR #支持向量机回归
from sklearn.neural_network import MLPRegressor #神经网络
#集成学习：Adaboost,随机森林,极端随机森林,提升树
from sklearn.ensemble import AdaBoostRegressor,RandomForestRegressor,ExtraTreesRegressor,GradientBoostingRegressor,BaggingRegressor
from xgboost import XGBRegressor #XGBoost

def regression_model(n_jobs):
    model_dict={}
    #线性回归
#     lr_model=LinearRegression(n_jobs=6) 
#     model_dict["线性回归"]=lr_model
    #岭回归
#     lars_model=Lars() 
#     model_dict["岭回归"]=lars_model
    #k近邻回归
    knr_model=KNeighborsRegressor(weights="uniform",n_jobs=13) 
    model_dict["k近邻回归"]=knr_model
    #回归树
    dtreer_model=DecisionTreeRegressor() 
    model_dict["回归树"]=dtreer_model
    #支持向量机回归
#     svr_model=SVR(kernel="linear") 
#     model_dict["支持向量机回归"]=svr_model
    #神经网络
    mlpr_model=MLPRegressor(max_iter=5000) 
    model_dict["神经网络"]=mlpr_model
    
    #Adaboost
    adaboost_model=AdaBoostRegressor() 
    model_dict["AdaBoost"]=adaboost_model
#     随机森林
    randomfr_model=RandomForestRegressor(n_jobs=n_jobs,random_state=10) 
    model_dict["随机森林"]=randomfr_model
    #极端随机森林
    etreer_model=ExtraTreesRegressor(n_jobs=n_jobs,random_state=10) 
    model_dict["极端随机森林"]=etreer_model
    #提升树
    gbr_model=GradientBoostingRegressor() 
    model_dict["提升树"]=gbr_model
    #XGBoost
    xgboostr_model=XGBRegressor(n_jobs=n_jobs,random_state=10) 
    model_dict["XGBoost"]=xgboostr_model
    #Bagging
#     bgr_model=BaggingRegressor(n_jobs=n_jobs,random_state=10) 
#     model_dict["BaggingRegressor"]=bgr_model

    return model_dict

"""分类模型"""
from sklearn.linear_model import LogisticRegression #Logistic
from sklearn.neighbors import KNeighborsClassifier #KNN
from sklearn.naive_bayes import MultinomialNB #朴素贝叶斯
from sklearn.tree import DecisionTreeClassifier #决策树
from sklearn.svm import SVC #支持向量机
from sklearn.neural_network import MLPClassifier #神经网络
#集成学习：Adaboost,随机森林,极端随机森林,提升树
from sklearn.ensemble import AdaBoostClassifier,RandomForestClassifier,ExtraTreesClassifier,GradientBoostingClassifier,BaggingClassifier
from xgboost import XGBClassifier #XGBoost
def classification_model(n_jobs):
    model_dict={}

#     lr_model=LogisticRegression(max_iter=1000) #Logistic
#     model_dict["Logistic"]=lr_model
    
#     knn_model=KNeighborsClassifier() #KNN
#     model_dict["KNN"]=knn_model
    
#     bys_model=MultinomialNB() #朴素贝叶斯
#     model_dict["朴素贝叶斯"]=bys_model
    
#     dtree_model=DecisionTreeClassifier() #决策树
#     model_dict["决策树"]=dtree_model
    
#     svm_model=SVC() #支持向量机
#     model_dict["支持向量机"]=svm_model
    
#     mlpc_model=MLPClassifier(max_iter=3000) #神经网络
#     model_dict["神经网络"]=mlpc_model

    adaboostc_model=AdaBoostClassifier(random_state=10)
    model_dict["AdaBoost"]=adaboostc_model

    randomfc_model=RandomForestClassifier(n_jobs=n_jobs,random_state=10)
    model_dict["随机森林"]=randomfc_model
    
    etreec_model=ExtraTreesClassifier(random_state=10) #极端随机森林
    model_dict["极端随机森林"]=etreec_model
    
    gbc_model=GradientBoostingClassifier(random_state=10) #提升树
    model_dict["提升树"]=gbc_model
    
    xgboostc_model=XGBClassifier(n_jobs=n_jobs,random_state=10) #XGBoost
    model_dict["XGBoost"]=xgboostc_model
    
    bgc_model=BaggingClassifier(n_jobs=n_jobs,random_state=10) #Bagging
    model_dict["Bagging"]=bgc_model
    
    return model_dict
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
