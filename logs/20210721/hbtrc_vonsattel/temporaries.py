import numpy as np
import torch

def findInteraction(featureList):
    #Generate linearized upper quadrant of interactions
    #Can generate lower quadrant by just modifying j's range from (i,n) to (i+1)
    #x = list(range(1,6))
    n = len(featureList)
    interaction = [None]*((n+1)*n//2)
    counter = 0
    for i in range(n):
        for j in range(i,n):
            interaction[counter] = featureList[i]*featureList[j]
            counter+=1
    #print(featureList)
    #print(interaction)
    return interaction

def findInteractionsSelf(featureMatrix):
    '''
    Generate linearized upper quadrant of interactions
    Can generate lower quadrant by just modifying j's range from (i,n) to (i+1)
    Input:
        featureMatrix: a torch tensor of shape [N,D] where N is number of samples and D is number of features
    Output:
        interaction: a torch tensor of shape [N,(D)(D+1)/2]
    '''
    n,d = featureMatrix.shape
    interaction = torch.zeros(n,((d+1)*d//2))
    counter = 0
    for i in range(d):
        for j in range(i,d):
            interaction[:,counter] = featureMatrix[:,i]*featureMatrix[:,j]
            counter+=1
    return interaction

def findInteractionsCross(feature1, feature2):
    '''
    Generate linearized matrix of interactions
    Input:
        feature1: a torch tensor of shape [N,D1] where N is number of samples and D1 is number of features of type 1
        feature2: a torch tensor of shape [N,D2] where N is number of samples and D2 is number of features of type 2
    Output:
        interaction: a torch tensor of shape [N, D1*D2]
    '''
    n,d1 = feature1.shape
    _,d2 = feature2.shape
    interaction = torch.zeros(n, d1*d2)
    counter = 0
    for i in range(d1):
        for j in range(d2):
            interaction[:,counter] = feature1[:,i]*feature2[:,j]
            counter+=1
    return interaction

'''
### Testing
x = np.array([[1,2,3,4],[5,6,7,8],[9,10,11,12]])
print(x.shape)
y = torch.FloatTensor(x)
print(y,y.shape)
print(findInteractions(y))
'''

##Train_test split
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit, train_test_split
def trte_split(filename):
    df = pd.read_csv(filename, delimiter=',')
    patient_id = df["patient_id"]
    indices = list(range(len(patient_id.values.tolist())))
    _,_,train_inds, test_inds = train_test_split(patient_id, indices, test_size=0.33, random_state=42)
    #print(set(patient_id[train_inds]).intersection(patient_id[test_inds]), patient_id[test_inds])
    #train_inds, test_inds = next(GroupShuffleSplit(test_size=.33, n_splits=2, random_state = 7).split(df, groups=df['cancer_subtype']))
    with open("BRCA_alt/trte_partition.txt", 'w') as file:
        file.write("training patient_id \n")
        file.write(",".join(patient_id[train_inds].values.tolist())+"\n")
        file.write("testing patient_id \n")
        file.write(",".join(patient_id[test_inds].values.tolist()))
    
    return train_inds, test_inds

#trte_split("BRCA_alt/PAM50_subtype.csv")


## Get label specific accuracies
def label_specific_acc(truth, prediction):
    """
    This function generates label specific accuracu scores.
    Input:
        truth:      a list of values constituting the ground truth
        prediction: a list of values constituting the predicted values by the model
    Output:
        accuracies: a dictionary with keys as the unique labels in the truth list and values as the label specific accuracies
    """
    assert len(truth) == len(prediction), "The the ground truth and prediction lists are not of the same length"

    labels = list(set(truth))
    num_labels = len(labels)
    total_counter = {label:0 for label in labels}
    correct_counter = {label:0 for label in labels}
    for i in range(len(truth)):
        total_counter[truth[i]] += 1
        if truth[i] == prediction[i]:
            correct_counter[truth[i]] += 1
    accuracies = {label:correct_counter[label]/total_counter[label] for label in labels}
    return accuracies

if __name__ == '__main__':
    ground_truth = [0,0,0,0,0,0,0,1,1,1,1,1,1,1,1] 
    prediction =  [0,0,0,0,0,0,0,0,1,1,1,1,1,1,0]
    label_specific_acc(ground_truth,prediction)