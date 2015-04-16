import sys
import os


__all__ = ["CNOfuzzy"]




class CNOfuzzy(object):
    """


    Optimise first time point in the cnolist and produces a report.




    - initial values are set to NA except for:
        - set initial values of the stimuli 
        - flip inhibitors to 1-inh (0 is inhibited and 1 is non-inhibited)
        - inhibitors where value is 1 are set to NA
        - initial values that are inhibited are set to the value of the inhibitors (1 where required)

      normHill_tf = function(x,g,n,k) {
          return(g*(1+k^n)*x^n/(x^n+k^n))
              }

       minNA<-function(x){
           if(all(is.na(x))){
               return(NA)
           }else{
               return(min(x,na.rm=TRUE))
           }
           }
       compOR<-function(x){
       if(all(is.na(x[which(simList$maxIx == s)]))){
       res<-NA
       }else{
       res<-max(x[which(simList$maxIx == s)],na.rm=TRUE)
       }
       return(res)
                      }


    For the simulation, set to NA dummies values if needed
    For each edge: 
        - compute normhill TF for each node
        - flip values that enter with a negative (1 - transval)

    Then, for each species, get the inputs and compute the compOR
    If an AND gate, compute the min. If all NA returns NA otherwise the min

    - reset inhibitors and stimuli. Here ceck original code !!
    - inhibitors = 1 - inhibitors

    - replace NAs with zeros to avoid having NA penalty applying to unconnected species
    
    after loop, set un-resolved bits to NA.


    """
    def __init__(self):
        pass
