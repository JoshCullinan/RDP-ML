import pandas as pd
import numpy as np
from pandas.core import series
from pandas.core.arrays import categorical
from pandas.core.indexes import category
from scipy.sparse import data
from sklearn import preprocessing
import sklearn 
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectFromModel
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay

#Input data
#acc = pd.read_csv(r'output/rdp_accuracy.txt', sep ='\t')

#def input(file_location):
file_location = 'output/ml_input.txt'
recom = pd.read_csv(file_location, sep='\t')
recom_phillip = pd.read_csv(r'output/ml_input_Phillip.txt', sep='\t')
recom = pd.concat([recom, recom_phillip])
#recom.drop(columns=['Consensus(A:0)', 'Consensus(A:1)', 'Consensus(A:2)'], inplace=True)

#Convert all entries from str to list
f = lambda x: x.strip('()')#.split(', ')
for (colname, coldata) in recom.iteritems():
    recom.loc(axis=1)[colname] = recom.loc(axis=1)[colname].apply(f)

Valid = recom.sample(frac=0.1).reset_index(drop=True)

#Convert data to single columns.
X = pd.DataFrame()
for i in recom.columns:
    current = pd.DataFrame(recom[i].str.split(',', expand=True).values, columns=[i, i, i])
    out = current.iloc[:,0].astype(float)
    out = out.append(current.iloc[:,1].astype(float), ignore_index=True)
    out = out.append(current.iloc[:,2].astype(float), ignore_index=True)
    X = X.join(out, how='outer')

X_valid = pd.DataFrame()
for i in Valid.columns:
    current = pd.DataFrame(Valid[i].str.split(',', expand=True).values)
    out = pd.Series(dtype=np.float64)
    for label, Series in current.iterrows():
        out = out.append(Series, ignore_index=True)
    out.name = i
    X_valid = X_valid.join(out, how='outer')


#Create Y
y = X.pop('Recombinant')
y_valid = X_valid.pop('Recombinant').astype(np.float64)

#Preprocessing
scaler = preprocessing.StandardScaler().fit(X)
X_scaled = pd.DataFrame(scaler.transform(X))
scaler = preprocessing.StandardScaler().fit(X_valid)
X_valid = pd.DataFrame(scaler.transform(X_valid))

#Make Test & Train -> Shuffles the test and train. Therefore, valid set made seperately. 
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, random_state=42)


#return X_train, X_test, y_train, y_test

#X_train, X_test, y_train, y_test = input(r'output/ml_input.txt')

acc = pd.read_csv('output/rdp_accuracy.txt',sep='\t')
ME = acc.loc(axis=1)['MatchedEvents'].sum()
RCC = acc.loc(axis=1)['RecombinantsCorrectlyChosen'].sum()
print(f'RDP ACCURACY TO BEAT: {(RCC/ME):.5f}')

### Model Fitting ###
#Logistic Regression
from sklearn.model_selection import learning_curve
clf = LogisticRegression(penalty = 'l2', 
                        random_state=42, 
                        multi_class='ovr', 
                        solver='liblinear', 
                        max_iter=1000)
clf.fit(X_train, y_train)
clf.score(X_test, y_test)

#Cross Validate Model
from sklearn.linear_model import LogisticRegressionCV
clf_cv = LogisticRegressionCV(Cs=10,
                            cv=5,
                            n_jobs=-1,
                            random_state=42,
                            penalty = 'l2', 
                            multi_class='ovr', 
                            solver='liblinear', 
                            max_iter=500)
clf_cv.fit(X_train, y_train)
clf_cv.score(X_test, y_test)

from sklearn.metrics import cohen_kappa_score
pred = clf_cv.predict(X_test)
cohen_kappa_score(y_test, pred)
disp = ConfusionMatrixDisplay(confusion_matrix(y_true=y_test, y_pred=pred), display_labels=['Parent', 'Recombinant'])
disp.plot()

#Learning Curve to visualise learning
import matplotlib.pyplot as plt
train_sizes, train_scores, test_scores, _, _ = learning_curve(LogisticRegression(), X_train, y_train, train_sizes=[0.1, 0.3,0.5,0.8,1.0], cv=5, return_times=True, n_jobs=-1)
train_scores_mean = np.mean(train_scores, axis=1)
train_scores_std = np.std(train_scores, axis=1)
test_scores_mean = np.mean(test_scores, axis=1)
test_scores_std = np.std(test_scores, axis=1)

plt.fill_between(train_sizes, train_scores_mean - train_scores_std,
                        train_scores_mean + train_scores_std, alpha=0.1,
                        color="r")
plt.fill_between(train_sizes, test_scores_mean - test_scores_std,
                        test_scores_mean + test_scores_std, alpha=0.1,
                        color="g")
plt.plot(train_sizes, train_scores_mean, 'o-', color="r",
                label="Training score")
plt.plot(train_sizes, test_scores_mean, 'o-', color="g",
                label="Cross-validation score")
plt.legend(loc="best")

#Test on test cases. Selecting the recombinant based on the strongest signal from the algorithm
pred_probs = clf_cv.predict_proba(X_valid)
predictions = pd.Series(dtype=np.float64)
for i in range(0, len(pred_probs)-1, 3):
    current = [pred_probs[i][1], pred_probs[i+1][1], pred_probs[i+2][1]]
    max_index = np.argmax(current)
    output = [0,0,0]
    output[max_index] = 1
    output = pd.Series(data=output)
    predictions = predictions.append(output, ignore_index=True)

disp = ConfusionMatrixDisplay(confusion_matrix(y_valid, predictions, display_labels=['Parent', 'Recombinant']))
disp.plot()

from sklearn.metrics import accuracy_score, classification_report
accuracy_score(y_valid,predictions)
print(classification_report(y_valid, predictions, target_names=['Parent', 'Recombinant']))

# #Feature selection from Model. -- Doesn't seem to help much.
# clf_featSelect = SelectFromModel(clf_cv, prefit=True)
# X_train_new = clf_featSelect.transform(X_train)
# X_test_new = clf_featSelect.transform(X_test)
# X_train_new.shape

# clf.fit(X_train_new, y_train)
# clf.score(X_test_new,y_test)


#Support Vector Machine
from sklearn import svm
clf2 = svm.SVC(max_iter=5000).fit(X_train, y_train)
clf2.score(X_test, y_test)

#Basic Neural Network
from sklearn.neural_network import MLPClassifier
clf3 = MLPClassifier(solver='adam', 
                    hidden_layer_sizes=(20,5), 
                    max_iter=35, 
                    random_state=42, 
                    #verbose=True, 
                    early_stopping=True, 
                    alpha = 0.01,
                    tol=0.0005,
                    n_iter_no_change= 15, 
                    activation='identity').fit(X_train, y_train)
clf3.score(X_test, y_test)


#Decision Trees
from sklearn.tree import DecisionTreeClassifier
clf4 = DecisionTreeClassifier(random_state=42)
clf4.fit(X_train, y_train)
clf4.score(X_test,y_test)


#Naive Bayes
from sklearn.naive_bayes import GaussianNB
clf5 = GaussianNB()
clf5.fit(X_train, y_train)
clf5.score(X_test, y_test)

#Extra Tress
from sklearn.ensemble import ExtraTreesClassifier
clf6 = ExtraTreesClassifier(n_estimators = 250, random_state=42)
clf6.fit(X_train, y_train)
clf6.score(X_test, y_test)

#Random Forest
from sklearn.ensemble import RandomForestClassifier
clf7 = RandomForestClassifier(n_estimators=250, random_state=42)
clf7.fit(X_train, y_train)
clf7.score(X_test,y_test)