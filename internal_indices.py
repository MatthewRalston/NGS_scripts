#!/usr/bin/env python
"""
------------------------------------------------------------------------------
--                                                                          --
--                                 MATTHEW RALSTON                          --
--                                                                          --
--                             S K E L E T O N . P Y                        --
--                                                                          --
------------------------------------------------------------------------------
-- TITLE                                                                    --
--                                                                          --
--  FallSpring 2014                                                         --
--                                                                          --
------------------------------------------------------------------------------
-- This file is designed to take a ordered set of features, an ordered list --
-- of their clusterings, and a distance metric. It calculates various       --
-- metrics for each cluster: the silhouette, davis-boudin, dunn coefficients--
-- the average dissimilarity, and min and max intercluster distance         --
--                                                                          --
------------------------------------------------------------------------------
"""
################################################
#
#               I M P O R T
#
################################################



import numpy as np
from scipy.spatial import distance as dis
import scipy as sp
import operator

################################################
#
#               U S E R    V A R I A B L E S
#
################################################



################################################
#
#               S U B - R O U T I N E S
#
################################################

################################################

################################################

def internal_indices(features,orderings,distance="euclidean"):
    clusters={}
    indices=[]
    for i,x in enumerate(orderings):
        if x in clusters.keys():
            clusters[x].append(features[i])
        else:
            clusters[x]=[features[i]]
    # 'A'
    centroids={}
    # 'B'
    avgdissim={}
    for i in clusters.keys():
        centroids[i]=np.mean(clusters[i],axis=0)
        sumdist=0
        for x in clusters[i]:
            sumdist+=eval("dis."+distance+"(x,centroids[i])")
        avgdissim[i]=sumdist/len(clusters[i])
    maxB=max(avgdissim)
    # 'D'
    dists={}
    for c,i in enumerate(clusters.keys()):
        # 'C'
        i_to_centroid=[]
        for j in np.delete(clusters.keys(),c):
            sumdist=0
            for x in clusters[i]:
                sumdist+=eval("dis."+distance+"(x,centroids[j])")
            i_to_centroid.append(sumdist/len(clusters[i]))
            a,b=sorted([i,j])
            dists[str(a)+str(b)]=eval("dis."+distance+"(centroids[i],centroids[j])")
        a,b=avgdissim[i],min(i_to_centroid)
        # average Silhouette of cluster
        silhouette=(b-a)/max([a,b])
        # Davies-Bouldin coefficient
        #    the average over all clusters would be the davies bouldin index
        temp=dict(avgdissim)
        del temp[i]
        d,e=max(temp.iteritems(),key=operator.itemgetter(1))
        a,b=sorted([i,d])
        dbc=(avgdissim[i]+e)/dists[str(a)+str(b)]
        # Dunn coefficient
        #     the minimum over all clusters would be the Dunn index
        temp=[s for s in dists.keys() if str(i) in s]
        temp={k: dists[k] for k in temp}
        di=min(temp.values())/maxB
        # add indices to the list
        a,b=min(temp.values()),max(temp.values())
        indices.append([silhouette,dbc,di,avgdissim[i],a,b])
    return indices

#*****************************************************************************#
################################################
#
#-----------------------------------------------
#
#                   M A I N
#-----------------------------------------------
#
################################################




##########################  E O F   ############################################
