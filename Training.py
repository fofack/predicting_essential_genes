
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split, cross_validate, RepeatedStratifiedKFold
from sklearn.ensemble import RandomForestClassifier
import joblib
from imblearn.ensemble import BalancedRandomForestClassifier
from statistics import mean
from imblearn.over_sampling import SMOTE
from sklearn import metrics
from Utils import Util
import pandas as pd
import csv
import os


class Training:

    def Training(dataset, group_name):
        X = dataset.drop('essentiality_status', axis=1)
        y = dataset['essentiality_status']

        # Use SMOTE to oversample the minority class
        oversample = SMOTE()
        over_X, over_y = oversample.fit_resample(X, y)
        over_X_train, over_X_test, over_y_train, over_y_test = train_test_split(
            over_X, over_y, test_size=0.2, stratify=over_y)
        # Build SMOTE SRF model
        SMOTE_SRF = RandomForestClassifier(n_estimators=150, random_state=0)
        # Create Stratified K-fold cross validation
        cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=1)
        scoring = ('f1', 'recall', 'precision')
        # Evaluate SMOTE SRF model
        scores = cross_validate(
            SMOTE_SRF, over_X, over_y, scoring=scoring, cv=cv)
        # #Get average evaluation metrics
        print('Mean f1: %.3f' % mean(scores['test_f1']))
        print('Mean recall: %.3f' % mean(scores['test_recall']))
        print('Mean precision: %.3f' % mean(scores['test_precision']))

        # Randomly spilt dataset to test and train set
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.1, stratify=y)
        # Train SMOTE SRF
        SMOTE_SRF.fit(over_X_train, over_y_train)
        # SMOTE SRF prediction result
        y_pred = SMOTE_SRF.predict(X_test)
        joblib.dump(value=SMOTE_SRF, filename='Models/'+group_name+'.pkl')
        results = classification_report(y_test, y_pred, output_dict=True)
        print(results)
        roc_auc = metrics.roc_auc_score(y_test, y_pred)
        print(roc_auc)

    def dataset_construction(input_path, subGroupName):
        data = pd.DataFrame()
        organism = Util.find_dir(Util(), input_path)
        for i in range(len(organism)):
            organismData = Util.find_dir(Util(), input_path+organism[i])
            for j in range(len(organismData)):
                if (organismData[j].split("_")[2].split(".")[0] == subGroupName):
                    tmp = pd.read_csv(
                        input_path+organism[i]+'/'+organismData[j])
                    data = pd.concat([data, tmp])
                    # print(len(data))
                    data.to_csv('data/'+subGroupName+'.csv', index=False)
        data = pd.DataFrame()
