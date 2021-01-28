# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 15:08:44 2019

@author: x
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 16:48:49 2017

@author: Lukasz 
"""

import itertools
import numpy as np
import matplotlib.pyplot as plt

from itertools import cycle
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_curve, auc
from scipy import interp

def PlotWithAnnotations(s, t, annIdx, annLabel, annColor, title, rpeaks = None, start=None, to=None):
    
    if(start is None):
        start = 0
    if(to is None):
        to = len(s)
    
    sCut = s[start:to]
    tCut = t[start:to]
    
    sMax = np.max(sCut)
    sMin = np.min(sCut)

    fig = plt.figure()
    plt.title(title)
    ax = fig.add_subplot(111)
    
    axes = fig.gca()
    axes.set_ylim([sMin, 1.5*max(sMax,0)])
    
    
    line, = ax.plot(tCut, sCut, lw=1)
        
    
        
    if rpeaks is not None:
        rpeaksIdx = t[rpeaks]
        rpeaksIdx = rpeaksIdx[np.logical_and( rpeaksIdx >= start, rpeaksIdx < to)]            
        ax.scatter( rpeaksIdx, np.ones_like(rpeaksIdx)*sMax*1.1, c='r')

    for i in range(len(annIdx)):
        if annIdx[i] >= start and annIdx[i] < to:
            xPos = tCut[annIdx[i]-start]
            ax.annotate(annLabel[i], xy = (xPos,1.2*sMax))  
            
            ax.vlines(x = xPos, ymax = 1.1*sMax, ymin = 0, color = annColor[i] )
    
    plt.show()   
    
    
def PlotConfusionMatrix(cm, classes,
                          normalize=False,
                          title='Confusion matrix',
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt),
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")

    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    
    
def PlotMulticlassRoc(y_test, y_score, class_names):
    n_classes = len(class_names)
    y_test_binarized = label_binarize(y_test, classes=class_names)

    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i, className in enumerate(class_names):
        fpr[className], tpr[className], _ = roc_curve(y_test == className, y_score[:, i])
        roc_auc[className] = auc(fpr[className], tpr[className])

        # Compute micro-average ROC curve and ROC area
    fpr["micro"], tpr["micro"], _ = roc_curve(y_test_binarized.ravel(), y_score.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr[i] for i in class_names]))

    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for className in class_names:
        mean_tpr += interp(all_fpr, fpr[className], tpr[className])

    # Finally average it and compute AUC
    mean_tpr /= n_classes

    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

    # Plot all ROC curves
    plt.figure()
    plt.plot(fpr["micro"], tpr["micro"],
         label='micro-average ROC curve (area = {0:0.2f})'
               ''.format(roc_auc["micro"]),
         color='deeppink', linestyle=':', linewidth=4)

    plt.plot(fpr["macro"], tpr["macro"],
         label='macro-average ROC curve (area = {0:0.2f})'
               ''.format(roc_auc["macro"]),
         color='navy', linestyle=':', linewidth=4)

    colors = cycle(['aqua', 'darkorange', 'cornflowerblue','green'])
    for i, color in zip(class_names, colors):
        plt.plot(fpr[i], tpr[i], color=color,
             label='ROC curve of class {0} (area = {1:0.2f})'
             ''.format(i, roc_auc[i]))

    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Some extension of Receiver operating characteristic to multi-class')
    plt.legend(loc="lower right")
    plt.show()
    return roc_auc
