# -*- python -*-
#
#  This file is part of cellnopt.core software
#
#  Copyright (c) 2011-2013 - EBI-EMBL
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: http://www.cellnopt.org
#
##############################################################################
import numpy
from scipy.optimize import fsolve
import pylab


__all__ = ["plot_ode_transfer_function", "normhill_tf",
"create_picture_cnorfuzzy_package"]



def plot_ode_transfer_function(tau, inputs):
    X = pylab.linspace(0,1,100)
    dX = ode_transfer_function(X, tau, inputs)
    from pylab import plot
    plot(X,dX)



def ode_transfer_function(x, tau, inputs):
    """

    :param list inputs: list of inputs as dictionaries that contains
        the ode parameters and nature of the link (activation or inhibition)

    .. plot::
        :include-source:
        :width: 80%

        from cellnopt.core.matfun import ode_transfer_function
        from pylab import linspace, plot, show
        X = linspace(0,2,100)
        dX = ode_transfer_function(X ,.1,
                {"B":{'k':1, 'n':.5, 'type':1},
                "Akt":{'k':0.5, 'n':2, 'type':-1}})
        plot(X,dX)
        show()

    """
    results = []
    for species, value in inputs.items():
        res = normhill_tf(x, value['k'], value['n'], g=1)
        if value['type'] == -1:
            res = 1 -res
        results.append(res)
    total = reduce(lambda x,y: x*y, results)
    dx = (total - x) * tau
    return dx


def normhill_tf(x, k, n, g=1):
    """Return the normalised Hill transfer function

    :param array x: Level of input data (from 0 to 1)
    :param float n: the Hill coefficient sharpness of the sigmoidal transition
        between high and low output node values
    :param float k: is the sensitivity parameter specyfying the EC50 value
    :param float g: a factor
    """
    return g * (k**n + 1) * x**n / (k**n + x**n)


def _hill_function(k, n=3, EC50=0.5):
    """
    inverse of the hill function
    from scipy.optimize import fsolve
    myres=fsolve(hill_function,0.5, args=(3,0.5))
    print('Results:'+str(myres))
    """
    return (1+k**n)*EC50**n/(EC50**n+k**n) -0.5


def _findEC50(x50,k,n):
    return x50**n * (1. + k**n) - 0.5*(x50**n + k**n)


def getEC50(k,n):
    myres = fsolve(_findEC50, .5, args=(k,n))
    return myres



def getk(N=[3,3,3,3,3,3,1.01], EC50=[0.2,0.3,0.4,0.5,0.6,0.7,0.5]):
    """How to get the k values if we know the EC50 and the n parameter"""
    k = []
    for n, ec50 in zip(N, EC50):
        myres = fsolve(_hill_function,0.5, args=(n, ec50))
        print('Results: n=%s' % n +str(myres))
        k.append(myres)
    return k


def create_picture_cnorfuzzy_package(fontsize=20, save=True):
    """Creates pictures of Fuzzy hill functions as in Morris et al.

    .. plot::
        :include-source:
        :width: 50%

        from cno.misc.matfun import create_picture_cnorfuzzy_package
        create_picture_cnorfuzzy_package()

    """
    n = [3, 3, 3, 3, 3, 3, 1.01]
    ec50 = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.5]
    k = getk(n, ec50)
    x = numpy.linspace(0, 1, 100)

    pylab.figure(1)
    pylab.clf()
    for i in range(0, len(n)):
        pylab.plot(x, normhill_tf(x,k[i],n[i]), linewidth=2,
                label="n=%s EC50=%s" % (n[i], ec50[i]))

    pylab.xlabel("Input Value", fontsize=22)
    pylab.ylabel("Output Value", fontsize=22)
    pylab.legend(loc="lower right")
    pylab.xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
            fontsize=fontsize)
    pylab.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
            fontsize=fontsize)

    if save:
        pylab.savefig('tf1.pdf')

    pylab.figure(2)
    pylab.clf()
    n=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    for i in range(0, len(n)):
        pylab.plot(x, x*n[i], linewidth=2)
    pylab.xlabel("Input Value", fontsize=22)
    pylab.ylabel("Output Value", fontsize=22)
    pylab.xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], fontsize=fontsize)
    pylab.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], fontsize=fontsize)
    if save:
        pylab.savefig('tf2.pdf')



def hill_activator(data, n=2, K=.5, betamax=1):
    r"""Hill function for an activator

    For an activator the Hill function is an increasing function with concentration of active activator X.
    The production rate rises from zero and reaches the maximal producion rate in a sigmoidal shape.

    .. math:: f(X) = (\beta_{max} X^n )/ (K^n + X^n)

    http://www.bio-physics.at/wiki/index.php?title=Hill_Function

    maximal production rate
    Hill coefficient
    activation coefficient

    :math:`\beta_{max}` is the maximal production rate of the promoter-transcirption factor complex.
    If the activator concentration X equals the activation coefficient K, the Hill function
    reaches the half-maximum point. The value of K is hardwired and a characteristic of the
    promoter, that depends on the promoter strength as well as on increased
    affinity through transcription factors.

    The lager the Hill coefficient n the more step like the Hill function becomes and in the
    limit n infinite the f(X) takes a unit step at X=K. The Hill coefficient comes from the fact
    the transcription factors can act as multimeres which leads to cooperative behaviour. Typical values for n are 1-4.

    """
    #if this is a dataframe, try:
    try:
        return betamax * (data**n).divide( data**n + np.array(K)**n)
    except:
        return betamax * (data**n / (data**n + K**n))

def hill_repressor(data, n=2, K=.5, betamax=1):
    r"""For repressors the Hill function decreases with the concentration of active repressor X

    .. math:: \beta_{max} / (1 + (X/K)**n)


    maximal production rate
    Hill coefficient
    activation coefficient

    The more repressor is available, the higher the probability that a repressor binds to the
    operator site, thus the expression level is more and more repressed with increasing
    repressor levels. Half-maximal repression occurs, when the concentration of active
    repressor equals the repression coefficient X=K


    """
    try:
        return betamax/(1. +  (data.divide(K))**n)
    except:
        return betamax/(1. +  (data/K)**n)



def normalise(df, n=2, Kact=.54, Krep=.26, tag_neg=True):
    fc = df.copy().astype(float)
    fc[df>0] = hill_activator(df, n=n, K=Kact)[df>0]
    if tag_neg is True:
        fc[df<0] = -1 * (1-hill_repressor(abs(df), n=n, K=Krep)[df<0])
    else:
        fc[df<0] = 1-hill_repressor(abs(df), n=n, K=Krep)[df<0]
    return fc






def hill_activator(data, n=2, K=.5, betamax=1):
    r"""Hill function for an activator

    For an activator the Hill function is an increasing function with concentration of active activator X.
    The production rate rises from zero and reaches the maximal producion rate in a sigmoidal shape.

    .. math:: f(X) = (\beta_{max} X^n )/ (K^n + X^n)

    http://www.bio-physics.at/wiki/index.php?title=Hill_Function

    maximal production rate
    Hill coefficient
    activation coefficient

    :math:`\beta_{max}` is the maximal production rate of the promoter-transcirption factor complex.
    If the activator concentration X equals the activation coefficient K, the Hill function
    reaches the half-maximum point. The value of K is hardwired and a characteristic of the
    promoter, that depends on the promoter strength as well as on increased
    affinity through transcription factors.

    The lager the Hill coefficient n the more step like the Hill function becomes and in the
    limit n infinite the f(X) takes a unit step at X=K. The Hill coefficient comes from the fact
    the transcription factors can act as multimeres which leads to cooperative behaviour. Typical values for n are 1-4.

    """
    #if this is a dataframe, try:
    try:
        return betamax * (data**n).divide( data**n + np.array(K)**n)
    except:
        return betamax * (data**n / (data**n + K**n))

def hill_repressor(data, n=2, K=.5, betamax=1):
    r"""For repressors the Hill function decreases with the concentration of active repressor X

    .. math:: \beta_{max} / (1 + (X/K)**n)


    maximal production rate
    Hill coefficient
    activation coefficient

    The more repressor is available, the higher the probability that a repressor binds to the
    operator site, thus the expression level is more and more repressed with increasing
    repressor levels. Half-maximal repression occurs, when the concentration of active
    repressor equals the repression coefficient X=K


    """
    try:
        return betamax/(1. +  (data.divide(K))**n)
    except:
        return betamax/(1. +  (data/K)**n)


def normalise(df, n=2, Kact=.54, Krep=.26, tag_neg=True):
    fc = df.copy().astype(float)
    fc[df>0] = hill_activator(df, n=n, K=Kact)[df>0]
    if tag_neg is True:
        fc[df<0] = -1 * (1-hill_repressor(abs(df), n=n, K=Krep)[df<0])
    else:
        fc[df<0] = 1-hill_repressor(abs(df), n=n, K=Krep)[df<0]
    return fc




