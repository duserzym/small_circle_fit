import pmag
from math import sqrt,atan,degrees,acos
import random

#============================================
#This points sort of make a cute little small circle:
#sclist=[ [220,43,1],[221,54,1],[229,37,1],[240,61,1],[243,36,1],[256,41,1],[259,55,1] ]

#Get eigenvectors of points at end of these unit vectors.

#"True" in arg list below means that we translate the center of mass of the vector endpoints
#(which are on the unit sphere) to the origin.

#Eigen finds the best great circle fit to the translated points. The normal to that great circle
#is also the center of the great circle translated out from the origin so the vector endpoints
#are back on the unit sphere.  The plane defined by the great circle (at the origin) cuts the unit sphere
#to form the small circle we are after.

#Lets see how it works for subsets of the 7 points

numpicks=int(input('How many of the 7 points do want to fit to? '))

trials=10
sclist=[]
for i in range(trials):
    biglist=[ [220,43,1],[221,54,1],[229,37,1],[240,61,1],[243,36,1],[256,41,1],[259,55,1] ]
    sclist=[]
    picklist=[]
    
    #randomly pick a subset, no dupes
    while len(picklist)<numpicks:
        pick=random.randint(0,6)
        if pick not in picklist:
            picklist.append(pick)

    for i in range(numpicks):       
        sclist.append(biglist[picklist[i]])
                      
    num=len(sclist)

    scmat=pmag.eigmat(sclist,True)    #load up the matrix
    eigvecs=pmag.eigen(scmat)         #get the eigenvectors.  We want the one with max eigenvalue.

    axdec=eigvecs[2][0]   
    axinc=eigvecs[2][1]
    axvec=[axdec,axinc,1]             #axis of the best fit small circle

    done=False


    #How far was the center of mass of vector endpoints from the origin?
    #The shorter the arc, the more this distance will larger than the distance
    #for a set of points forming a longer arc on the same great circle.
    nsum=0
    esum=0
    dsum=0
    for i in range(num):
        n,e,d=pmag.nedlist(sclist[i])
        nsum=nsum+n
        esum=esum+e
        dsum=dsum+d

    nav=nsum/num
    eav=esum/num
    dav=dsum/num

    dist=sqrt(nav*nav+eav*eav+dav*dav)
    #print('dist',dist)

    #because dist is too short for small arcs, this radius will be too low 
    radius=degrees(acos(dist/1.0))
    #print('radius',radius)    

    #Probably better to just do it this way:
    #Get list of angles from point to small circle axis, av angle, and av of squared angles
    anglist=[]
    angsum=0
    angsumsq=0
    N=len(sclist)

    for i in range(num):
        ang=pmag.ang(axdec,axinc,sclist[i][0],sclist[i][1])
        anglist.append(ang)
        angsum=angsum+ang
        angsumsq=angsumsq+ang*ang
        #print (ang)

    avang=angsum/num
    avang2=sqrt(angsumsq/num)

    #print('Three estimates of radius. For short arcs, first will be low.  {0:5.1f} {1:5.1f} {2:5.1f}'.format(radius,avang,avang2)) 
    #print('Probably best to use last one - sqrt of av of squared angles between data and circle center')
    print()
    print(picklist)
    print ('Best fit small circle:  Dec={0:5.1f} Inc={1:5.1f} Radius={2:5.1f}'.format(axdec,axinc,avang2))
    #print ('Av of angles (or angles-sq) from small circle center to data:  {0:5.1f} {1:5.1f}'.format(avang,avang2))


    
