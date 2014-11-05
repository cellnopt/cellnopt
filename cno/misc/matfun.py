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
from pylab import linspace

__all__ = ["ode_transfer_function", "normhill_tf",
"create_picture_cnorfuzzy_package"]



def plot_ode_transfer_function(tau, inputs):
    X = linspace(0,1,100)
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
    for species, value in inputs.iteritems():
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


def create_picture_cnorfuzzy_package():
    """Creates pictures of Fuzzy hill functions as in Morris et al.


    .. plot:: 
        :include-source:
        :width: 50%

        from cellnopt.core.matfun import create_picture_cnorfuzzy_package
        create_picture_cnorfuzzy_package()

    """
    from pylab import plot, legend, xlabel, ylabel, xticks, yticks,clf, savefig
    n=[3,3,3,3,3,3,1.01]
    ec50=[0.2,0.3,0.4,0.5,0.6,0.7,0.5]
    k = getk(n, ec50)
    x = numpy.linspace(0,1,100)
    clf()
    for i in range(0, len(n)):
        plot(x, normhill_tf(x,k[i],n[i]), linewidth=2, label="n=%s EC50=%s" % (n[i], ec50[i]))
    xlabel("Input Value", fontsize=22)
    ylabel("Output Value", fontsize=22)
    legend(loc="lower right")
    xticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1], fontsize=20)
    yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],fontsize=20)
       
    savefig('tf1.pdf')


    clf()
    n=[0.2,0.3,0.4,0.5,0.6,0.7,0.8]
    for i in range(0, len(n)):
        plot(x, x*n[i], linewidth=2)
    xlabel("Input Value", fontsize=22)
    ylabel("Output Value", fontsize=22)
    xticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1], fontsize=20)
    yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],fontsize=20)
       
    savefig('tf2.pdf')




