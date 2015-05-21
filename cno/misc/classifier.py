
import os
import numpy as np
import pandas as pd
from cno import *
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
import pylab


class Classifier(object):
    def __init__(self, midas=None):
        self.midas = midas

    def create_class(self, ds, thres, classes4=False):
         m0 = np.logical_and(ds.query('time==0') == 0. , ds.query('time!=0') > thres)
         ds[m0] = 0
 
         m1 = np.logical_and(ds.query('time==0') == 1. , ds.query('time!=0') <= thres)
         ds[m1] = 1
 
         m2 = np.logical_and(ds.query('time==0') == 0. , ds.query('time!=0') <= thres)
         ds[m2] = 2
         if classes4 is True:
             m3 = np.logical_and(ds.query('time==0') == 1. , ds.query('time!=0') > thres)
             ds[m3] = 3
         else:
             m3 = np.logical_and(ds.query('time==0') == 1. , ds.query('time!=0') > thres)
             ds[m3] = 2
 
         classes = ds.query('time==0')
         return classes

    def compute_stats(self, thres=0.5, classes4=False):
        confusion = {}
        classification = {}
        acc = {}
        prec_pos = {}
        prec_none = {}
        rec_pos = {}
        rec_none = {}
        f_pos = {}
        f_none = {}

        c_df = self.midas.df.copy()
        c_sim = self.midas.sim.copy()
        real = self.create_class(c_df, thres, classes4)
        pred = self.create_class(c_sim, thres, classes4)

        self.real = real
        self.pred = pred

        if classes4 is True:
            target_names = ['class 0', 'class 1', 'class 2', 'class 3']
            labels = [0, 1, 2, 3]
        else:
            target_names = ['class 0', 'class 1', 'class 2']
            labels = [0, 1, 2]

        for spp in c_df.columns:

            x = real[spp].values
            y = pred[spp].values
            m1 = np.isnan(x) == False
            m2 = np.isnan(y) == False
            mask = np.logical_and(m1, m2)
            x = x[mask]
            y = y[mask]

            classification[spp] = classification_report(x, y,
                                                            labels=labels,
                                                            target_names=target_names)
            confusion[spp] = confusion_matrix(x,y, labels=labels)
            acc[spp] = accuracy_score(x.flatten(), y.flatten())
            # Precision, recall and F1 calculated for positive FC (*_pos) and None FC (*_none)
            xred = x
            xred[xred == 3] = 2
            yred = y
            yred[yred==3] = 2

            prec_pos[spp] = precision_score(xred, yred, labels=[0, 1, 2], average=None)[0]
            prec_none[spp] = precision_score(xred, yred, labels=[0, 1, 2], average=None)[2]

            rec_pos[spp] = recall_score(xred, yred, labels=[0, 1, 2], average=None)[0]
            rec_none[spp] = recall_score(xred, yred, labels=[0, 1, 2], average=None)[2]

            f_pos[spp] = f1_score(xred, yred, labels=[0, 1, 2], average=None)[0]
            f_none[spp] = f1_score(xred, yred, labels=[0, 1, 2], average=None)[2]

        self.classification = classification
        self.accuracy = acc
        self.confusion = confusion
        self.precision_pos = prec_pos
        self.precision_none = prec_none
        self.recall_pos = rec_pos
        self.recall_none = rec_none
        self.f1_pos = f_pos
        self.f1_none = f_none

        # for all as well
        x = real.values
        y = pred.values
        m1 = np.isnan(x) == False
        m2 = np.isnan(y) == False
        mask = np.logical_and(m1, m2)
        x = x[mask]
        y = y[mask]

        self.classification_all = classification_report(x.flatten(), y.flatten(),
                                               labels=labels, target_names=target_names)
        self.confusion_all = confusion_matrix(x.flatten(), y.flatten(), labels=labels)
        self.accuracy_all = accuracy_score(x.flatten(), y.flatten())

        xred = x
        xred[xred == 3] = 2
        yred = y
        yred[yred==3] = 2

        self.precision_all_pos = precision_score(xred, yred, labels=[0, 1, 2], average=None)[0]
        self.precision_all_none = precision_score(xred, yred, labels=[0, 1, 2], average=None)[2]

        self.recall_all_pos = recall_score(xred, yred, labels=[0, 1, 2], average=None)[0]
        self.recall_all_none = recall_score(xred, yred, labels=[0, 1, 2], average=None)[2]

        self.f1_all_pos = f1_score(xred, yred, labels=[0, 1, 2], average=None)[0]
        self.f1_all_none = f1_score(xred, yred, labels=[0, 1, 2], average=None)[2]

        if classes4 is True:
            TP_pos = self.confusion_all[0, 0]
            FP_pos = np.sum(self.confusion_all[1:4, 0])
            P_pos = np.sum(self.confusion_all[0])
            N_pos = np.sum(self.confusion_all[1:4])
            FN_pos = np.sum(self.confusion_all[0, 1:4])

            TP_0 = np.sum(self.confusion_all[2:4, 2:4])
            FP_0 = np.sum(self.confusion_all[0:2, 2:4])
            P_0 = np.sum(self.confusion_all[2:4])
            N_0 = np.sum(self.confusion_all[0:2])
            FN_0 = np.sum(self.confusion_all[2, 0:2])

        else:
            TP_pos = self.confusion_all[0, 0]
            FP_pos = np.sum(self.confusion_all[1:3, 0])
            P_pos = np.sum(self.confusion_all[0])
            N_pos = np.sum(self.confusion_all[1:3])
            FN_pos = np.sum(self.confusion_all[0, 1:3])

            TP_0 = self.confusion_all[2, 2]
            FP_0 = np.sum(self.confusion_all[0:2, 2])
            P_0 = np.sum(self.confusion_all[2])
            N_0 = np.sum(self.confusion_all[0:2])
            FN_0 = np.sum(self.confusion_all[2:4, 0:2])

        self.TPR_pos = TP_pos/float(P_pos)
        self.FPR_pos = FP_pos/float(N_pos)

        self.TPR_none = TP_0/float(P_0)
        self.FPR_none = FP_0/float(N_0)


    def plot_confusion(self, species=None, cmap=None, tight_layout=False):
        if cmap is None:
            import colormap
            cmap = colormap.cmap_builder('white',
                    'blue','darkblue')
        from biokit import imshow
        if species is not None:
            imshow(self.confusion[species], cmap=cmap)
            pylab.title(species,fontsize=20)
            if tight_layout is True:
                pylab.tight_layout()

        else:
            imshow(self.confusion_all, cmap=cmap)


