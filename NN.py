import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from sklearn import preprocessing
from sklearn import metrics
from sklearn.model_selection import train_test_split
from sklearn.utils.validation import column_or_1d
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader

torch.cuda.is_available()

#Hyperparameters
EPOCHS = 80
BATCH_SIZE = 32
LEARNING_RATE = 0.001


file_location = 'output/ml_input.txt'

recom = pd.read_csv(file_location, sep='\t')
recom_phillip = pd.read_csv(r'output/ml_input_Phillip.txt', sep='\t')
recom = pd.concat([recom, recom_phillip])
#recom.drop(columns=['Consensus(A:0)', 'Consensus(A:1)', 'Consensus(A:2)'], inplace=True)

acc = pd.read_csv('output/rdp_accuracy.txt',sep='\t')
ME = acc.loc(axis=1)['MatchedEvents'].sum()
RCC = acc.loc(axis=1)['RecombinantsCorrectlyChosen'].sum()
print(f'RDP ACCURACY TO BEAT: {(RCC/ME):.5f}') 

#Convert all entries from str to list
f = lambda x: x.strip('()')#.split(', ')
for (colname, coldata) in recom.iteritems():
    recom.loc(axis=1)[colname] = recom.loc(axis=1)[colname].apply(f)

#Convert data to single columns.
X = pd.DataFrame()
for i in recom.columns:
    current = pd.DataFrame(recom[i].str.split(',', expand=True).values, columns=[i, i, i])
    out = current.iloc[:,0].astype(float)
    out = out.append(current.iloc[:,1].astype(float), ignore_index=True)
    out = out.append(current.iloc[:,2].astype(float), ignore_index=True)
    X = X.join(out, how='outer')

#Create Y
y = X.pop('Recombinant')
y = y.astype('category')
encode_map = {
    1.0: 1,
    0.0: 0
}
y.replace(encode_map, inplace=True)


#Preprocessing
scaler = preprocessing.StandardScaler().fit(X)
X_scaled = scaler.transform(X)
X_scaled.shape
#Make Test & Train
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, random_state=42, test_size=0.1, shuffle=True)

#Feature selection
from sklearn.linear_model import LogisticRegressionCV
from sklearn.feature_selection import SelectFromModel
from sklearn.feature_selection import RFE
clf_cv = LogisticRegressionCV(Cs=10,
                            cv=5,
                            n_jobs=-1,
                            random_state=42,
                            penalty = 'l2', 
                            multi_class='ovr', 
                            solver='liblinear', 
                            max_iter=100)
clf_cv.fit(X_train, y_train)
clf_cv.score(X_test, y_test)

#Feature selection from Model.
clf_featSelect = SelectFromModel(clf_cv, prefit=True)
clf_featSelect2 = RFE(clf_cv, n_features_to_select= 0.5, verbose = 1)
clf_featSelect2= clf_featSelect2.fit(X_train, y_train)
X_train = clf_featSelect2.transform(X_train)
X_test = clf_featSelect2.transform(X_test)
X_train.shape

clf_cv.fit(X_train, y_train)
clf_cv.score(X_test,y_test)


#DataLoader
## train data
class trainData(Dataset):
    
    def __init__(self, X_data, y_data):
        self.X_data = X_data
        self.y_data = y_data
        
    def __getitem__(self, index):
        return self.X_data[index], self.y_data[index]
        
    def __len__ (self):
        return len(self.X_data)
train_data = trainData(torch.FloatTensor(X_train), 
                       torch.FloatTensor(y_train))

#train_data = trainData(X_train, y_train)
## test data    
class testData(Dataset):
    
    def __init__(self, X_data):
        self.X_data = X_data
        
    def __getitem__(self, index):
        return self.X_data[index]
        
    def __len__ (self):
        return len(self.X_data)
    
test_data = testData(torch.FloatTensor(X_test))

train_loader = DataLoader(dataset=train_data, batch_size=BATCH_SIZE, shuffle=True)
test_loader = DataLoader(dataset=test_data, batch_size=3)

### CLASS METHOD ###
class binaryClassification(nn.Module):
    def __init__(self):
        super(binaryClassification, self).__init__()
        # Number of input features is 12.
        self.layer_1 = nn.Linear(20, 128) 
        self.layer_2 = nn.Linear(128, 64)
        self.layer_out = nn.Linear(64, 1) 
        
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(p=0.1)
        self.batchnorm1 = nn.BatchNorm1d(128)
        self.batchnorm2 = nn.BatchNorm1d(64)
        
    def forward(self, inputs):
        x = self.relu(self.layer_1(inputs))
        x = self.batchnorm1(x)
        x = self.relu(self.layer_2(x))
        x = self.batchnorm2(x)
        x = self.dropout(x)
        x = self.layer_out(x)
        
        return x

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

model = binaryClassification()
model.to(device)
print(model)
criterion = nn.BCEWithLogitsLoss()
optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE)


def binary_acc(y_pred, y_test):
    y_pred_tag = torch.round(torch.sigmoid(y_pred))

    correct_results_sum = (y_pred_tag == y_test).sum().float()
    acc = correct_results_sum/y_test.shape[0]
    acc = torch.round(acc * 100)
    
    return acc

model.train()
for e in range(1, EPOCHS+1):
    epoch_loss = 0
    epoch_acc = 0
    for X_batch, y_batch in train_loader:
        X_batch, y_batch = X_batch.to(device), y_batch.to(device)
        optimizer.zero_grad()
        
        y_pred = model(X_batch)
        
        loss = criterion(y_pred, y_batch.unsqueeze(1))
        acc = binary_acc(y_pred, y_batch.unsqueeze(1))
        
        loss.backward()
        optimizer.step()
        
        epoch_loss += loss.item()
        epoch_acc += acc.item()
        

    print(f'Epoch {e+0:03}: | Loss: {epoch_loss/len(train_loader):.5f} | Acc: {epoch_acc/len(train_loader):.3f}')

y_pred_list = []
model.eval()
with torch.no_grad():
    for X_batch in test_loader:
        X_batch = X_batch.to(device)
        y_test_pred = model(X_batch)
        y_test_pred = torch.sigmoid(y_test_pred)
        y_pred_tag = torch.round(y_test_pred)
        y_pred_list.append(y_pred_tag.cpu().numpy())

y_pred_list = [a.squeeze().tolist() for a in y_pred_list]

from sklearn.metrics import confusion_matrix, classification_report
confusion_matrix(y_test, y_pred_list)
print(classification_report(y_test, y_pred_list))


y_pred_list = []
model.eval()
with torch.no_grad():
    for X_batch in test_loader:
        X_batch = X_batch.to(device)
        y_test_pred = model(X_batch)
        y_test_pred = torch.sigmoid(y_test_pred)
        y_pred_tag = torch.round(y_test_pred)
        y_pred_list.append(y_pred_tag.cpu().numpy())

y_pred_list = [a.squeeze().tolist() for a in y_pred_list]
y_pred_list
y_preds = []
y_preds.append(y_test_pred.cpu().numpy())
y_preds = [a.squeeze().tolist() for a in y_preds]

X_test = X_test[:2667,:]
y_test = y_test[:2667]

out = []
model.eval()
with torch.no_grad():
    for X_batch in test_loader:
        X_batch = X_batch.to(device)
        y_test_pred = model(X_batch)
        y_test_pred = torch.sigmoid(y_test_pred).cpu().numpy()
        max = -1*math.inf
        for i in y_test_pred:
            if i > max:
                max = i
       
        for i in y_test_pred:
            if i == max:
                out.append(1)
            else:
                out.append(0)
out
test = pd.DataFrame(columns=[0,1,2])
while len(out) > 0:
    A = out.pop(0)
    A
    B = out.pop(0)
    C = out.pop(0)
    test.append([A,B,C])


confusion_matrix(y_test, out)
print(classification_report(y_test, out))