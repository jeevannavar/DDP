import seaborn as sns
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

def conf_matrix(y_label, y_pred, features=None):
    """
    y_label: True labels
    y_pred: Predicted labels
    features: labels along x-axis/y-axis of the confusion matrix
    plot confusion matrix
    """
    assert len(y_pred) == len(y_label)
    if features:
        assert len(np.unique(y_label)) == len(features)
    num_classes = max(len(np.unique(y_pred)), len(np.unique(y_label)))
    mat = np.array([[0 for i in range(num_classes)] for j in range(num_classes)])
    for i in range(len(y_pred)):
        mat[round(y_pred[i])][round(y_label[i])] +=1
    
    output_class_acc = [round(mat[i][i]/np.sum(mat[i]) *100,2) for i in range(num_classes)];
    target_class_acc = [round(mat[i][i]/np.sum(mat[:,i]) *100,2) for i in range(num_classes)]

    percentages = [[(mat[i][j]/len(y_pred))*100 for j in range(num_classes)] for i in range(num_classes)]
    labels = [[str(mat[i][j])+"\n"+str(round(percentages[i][j],2))+"%" for i in range(num_classes)] for j in range(num_classes)]
    
    percentages.append(target_class_acc)
    labels.append([str(output_class_acc[i]) + "%\n" + str(100 - output_class_acc[i]) for i in range(num_classes)])
    for i in range(num_classes):
        percentages[i].append(output_class_acc[i])
        labels[i].append(str(target_class_acc[i]) + "%\n" + str(100 - target_class_acc[i]))
    total = 0;
    for i in range(num_classes):
        total += 100*mat[i][i]/len(y_pred)
    percentages[-1].append(total)
    total = round(total, 2)
    labels[-1].append(str(total) + "%\n" + str(100 - total))

    colors = [plt.cm.tab20(0),plt.cm.tab20(1),plt.cm.tab20c(4),
          plt.cm.tab20c(5),plt.cm.tab20c(6),plt.cm.tab20c(7)]
    cmap=matplotlib.colors.ListedColormap(colors)

    sns.heatmap(percentages, vmin = 0, vmax = 100, annot = False, cmap = cmap)
    plt.xlabel("Target Class")
    plt.ylabel("Output Class")
    plt.title("Confusion Matrix")
    if features:
        plt.xticks(np.arange(0,len(features))+0.5, features + [''],rotation = 90)
        plt.yticks(np.arange(0,len(features))+0.5, features + [''], rotation = 0)

    for i in range(num_classes+1):
        for j in range(num_classes+1):
            text = plt.text(j+0.5, i+0.5, labels[j][i], ha="center", va="center", color="white", fontsize = 10, fontweight = 'bold')
    plt.show()

if __name__ == "__main__":
    l = [0 for i in range(20)]+ [1 for i in range(50)] + [2 for i in range(30)] # label
    m = [0 for i in range(43)] + [1 for i in range(10)] + [2 for i in range(47)] # prediction
    conf_matrix(l,m)