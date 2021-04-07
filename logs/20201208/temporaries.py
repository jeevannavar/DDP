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

