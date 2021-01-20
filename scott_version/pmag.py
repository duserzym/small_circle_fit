#===========================================================
# pmag module   
# BluntObject Software
# "It ain't pretty, but it usually runs"
#
# Useful fuctions for Paleomagnetism Practitioners
# who are writing Assorted Small Programs (ASPs)
# or Big Ol'Applications (BOAs) in Python.
#
#Ported to python 3.7 4 Oct 2020
#
#Bug(?)in getDataFromRMG rooted out 1/9/16
#6 Jan 2017  added orthoreg and allreg for line slopes
#22 Mar 2019 added coord conversions local to geo, geo to local
#Jan 2020 eigen2x2
#Jan 2020 matrix stuff used in getting AMS code going

#ALSO:  see PYTHON CHEAT SHEET BELOW!
#
#I/O stuff
#  getdirlist()  Creates a text file of file names in selected dir
#  GetFileNames(SamFileName) returns [Sam file name, path, list of files]
#  ListFromTextFile()   returns list of text lines
#  ListFromTextFile2()  returns [line,path,filename] so you can cwdir to folder
#  ListFromTextFileExt(ext)  same, but filter on ext like "rmg" or "txt"
#  ListToTextFile(list of text lines)
#  getNums(Textline,n)  returns list of n nums
#  getAllItems(comma-delimited-text-line)   returns list of nums (text if error)
#  NumListToText([list of nums])  returns line of text + CR
#  getDIMA[textline from Oxy data file]  returns ['D','I','M','ang','tag']
#  getDIMAnums[textline from Oxy data file]  returns [D, I, M, ang, 'tag']
#  getstratDIMAnums   same as above, but the tilt-corrected data
#  getDataFromRMG(biglist)  returns list of lists [tag,treat,Mz,D,I,M,ARM bias,secs since midnight]
#  
#Stats, mostly Fisher
#  FishStat(N,R)   returns [k,a95]
#  FakeFish(trueD, trueI, k, N)   returns list of N [D,I] pairs
#  Fisher(list of [D,I]pairs)   returns [meanD, meanI, k, a95,R]
#  mixfit(dirlist,polelist,guess)  dir+demag circles. Returns [D,I,R,k,a95,N,iter,angsum)
#  delvgp(list of [lon,lat] pairs)   returns delvgp
#  getasd (list[lon,lat] pairs)  returns [ASD-mean, ASD-spinaxis, mean lon, mean lat]
#  DecIncStats(list of [D,I] pairs)  returns [meanD, dsd, MeanI, isd]
#  booter(list)  returns resampled list
#  
#GP stuff
#  PaleoDir(vgp-lon, vgp-lat,site-lon, site-lat)  returns [D,I]
#  vgp(D,I,site lon, site lat) returns [vgp-lon, vgp-lat]
#  vgperr(D,I,slon,slat,a95) returns [vlon,vlat,dp,dm,b95]
#  kmodG(lat)    returns k   for 0-5 Ma
#  kmodGmio(lat) returns k   for Miocene 5-22.5 Ma
#  Flatten (I,f) returns I-flat
#  UnFlatten (I,f)  returns I-unflat
#  GetField(Lon,Lat,model)returns (dec,inc,H,x,y,z)  needs coeff. file
#      model params in WMM format.  Pass model file as model.    
#  GetField2010 (lon,lat,model) returns (dec,inc,H,x,y,z)  needs coeff. file
#      model params --  1=WMM "WMM.COF"  2=NAD "WMM-NAD.COF" 3=NDF  "WMM-NDF.COF"
#
#Rotation
#  HRot(D,I,azm,rot)  returns [D,I]  rot is CW
#  ARot(D,I,azm,plg,rot)  returns [D,I]  rot is CW
#  doyz[D,I,y,z]  returns [D,I] corrected for y,z  (strike, dip of core top)
#  DoYZ  same thing
#  undoyz[D,I,y,z]  returns [D,I] from geogr back to spec coords
#
#Vectors
#  ned(D,I,M)  returns [N,E,D]
#  dim(N,E,D)  returns [D,I,M]
#  nedlist([D,I,M]) list input; returns [N,E,D]
#  dimlist([N,E,D]) list input; returns [D,I,M]
#  veclen([N,E,D])  returns length
#  dot([N,E,D],[N,E,D]) returns [dotnum, ang]
#  ang(Dec1,Inc1,Dec2,Inc2) returns angle between them in degrees
#  cross ([N,E,D],[N,E,D]) returns [N,E,D] 
#  crossunit ([N,E,D],[N,E,D]) returns [N,E,D] as unit vector
#  VecSum(list of [D,I,M])  returns D,I,R
#  decomp[D,I,M,D1,I1,D2,I2]  returns [mom1,mom2,projq,DecompGood]
#  localtogeo[sitelon,sitelat,dec,inc]
#  geotolocal[sitelon,sitelat,lon,lat(describing the geocentric vector)]
#
#Matrix stuff
#  matmaker(list)   Makes an array (list of row-lists) from a list
#  matunmaker(array)  Takes an a array and returns a list
#  getrotmat(axisdec,axisinc,rot)   Right-hand, CCW rotation
#  dorotmat(dec,inc,axisdec, axisinc, rot)  returns [dec,inc,n,e,d]
#  ----------- General purpose matrix routines from web-------
#  tranposeMatrix(matrix)
#  getMatrixMinor(matrix,i,j)
#  getMatrixDeterminant(matrix)
#  getMatrixInverse(matrix)
#  matmult(matrix1,matrix2)    matrix1 X matrix2
#  ------------------------------------------------------------
#
#Line and plane fits  etc
#  eigen([k11,k12,k13,k22,k23,k33])  returns [[D,I,eigval] x 3],ok]
#  eigmat(list of [D,I,M] sets)  returns [nn,ne,nd,ee,ed,dd]
#  eigen2x2 Wants [[11,12],[21,22]] Gives [[11,12,ev1],[21,22,ev2]]  unit vecs
#  linefit(list of [D,I,M,A], free)  returns [fitD, fitI, mad, pathindex]
#  planefit(list of [D,I])   return [poledec, poleinc,eigenvalue, n]
#  closest(dec,inc,poledec,poleinc)  gives pt on gc closest to dir [gcdec,gcinc,ang]
#  orthoreg(list of x vals, list of y vals)  returns (slope,y-intercept)
#  allreg(list of x vals, list of y vals)  returns (m,my,mx))   3 slopes
#  ListStats(list of vals)  returns [mean,sd, sd/mean,n]
#  thellfit(xvals,yvals)  returns [b,sb,sb/b,N,Yint]  slope and its error
#
#Miscellaneous
#  GetLev(min,max,n) returns list of n (odd!) contour levels, half>0
#
#==========================================================
#  PYTHON CHEAT SHEET  8/18
#  --------------------
#  See DISLIN hints below
#
#  To continue a line, use backslash \
#  or put the break inside a (),[], or {} statment (preferred)
#
#  To see default paths for python to find files:
#  import sys
#  sys.path
#
#  To save path...
#       path=os.path.dirname(a)
#  ... and start there next time.
#       os.chdir(path)

#
#  To make it so that .py file can run by just clicking the file:
#  import os
#  os.system("pause")   put this at end of program to keep window open
#
#  For help in IDLE, type >>> help (file)
#
#  Generic while loop--
#  i=0
#  while i<num:
#     "do something"
#     i=i+1
#
#  Generic for loop (i will be 0, 1, 2, 3) --
#  for i in range(4):
#     do something
#
#  In loops:
#     break - takes you out of loop
#     continue - short circuits loop (but continues interating)
#     pass - placeholder after colon
#
#  If you don't need the index, for loops over items automtically!
#  for DummyItemName in (mylist):
#      print DummyItemName   #will print items in order!?
#
#  If you want the index,pythonically:
#  for i,item in enumerate(mylist):    #assign both i and item at once!
#      print i
#      print item
#
#  To just catch the enumerated list,  a=list(enumerate(mylist))
#
#  Generic if/else statment (Note True and False are predefined)
#    MyBool= True
#    if MyBool:
#      "do something
#    else:
#      "do this instead"
#
#  range(5) means [0,1,2,3,4]
#
#  To make a list: thisList= [one,two,three]  note: square brackets
#  ---  With parentheses, it's a tuple  (1,2,3)
#
#  Make an empty list to append to:    NewList=list()
#
#  Append new item to list:    NewList.append(newItem)
#  To get number of lines:  NumLines=len(ListOfLines)
#
#  To make a list via list comprehension
#   s=[2*x for x in range(5)]  gives [0,2,4,6,8]
#   s=[x*x for x in range(10) if 3*x>15]   gives [36,49,64,81]
#
#  Unpack a tuple or list:  a,b,c,d=[1,2,3,4]
#
#  a=[1,2,3]
#  a=b  CAREFUL!  a and b are just 2 names ('aliases') for the same list
#  Change a, b is changed too!
#
#  To clone a list   a=b[:]     (whole slice of list)
#  Now these are two different lists.  Watch out if it's a list of lists!
#
#  FUN WITH STRINGS (OR LISTS)
#  word[:7] slices out first 7 chars (index 0 up to but not incl. index 7)
#  word[7:] slices out index 7 to end
#  word[2:8] gives you 6 chars (from index 2 up to but not incl. 8)
#  word[:-1] gives all but last 
#
#  Formatted print: print('a=%.1f b=%5.1f c=%15.5f' %(1,2,3))
#      (would give      a=1.0 b=  2.0 c=       3.00000)
#
#  Or, use the new string method .format  Note - uses curly brackets:
#    print "{0:6.4f}  {1:.4E} {2:4} {3:04}".format (3.141456, 79.8, 4, 4)
#      3.1415  7.9800E+01    4 0004
#
#  To get incs (say) printing in 4 spaces, either -34.9  or [ ]34.9, use -5.1f  
#
#  To get 3 decimal places with no leading spaces:  0.3f  '3.457'and '-159.457'
#
#  More here:  http://docs.python.org/library/string.html#formatstrings
#
#  To get text from file into list: ListOfLines=OpenedFile.readlines()
#
#  To modify a filename or ext:
#    fname = os.path.basename(a)
#    part1=os.path.splitext(fname)[0]+'WHATEVER YOU MIGHT WANT TO ADD'
#    ext=os.path.splitext(fname)[1]
#    fname=part1+ext
#    dirname=os.path.dirname(a)
#    fullpathname=dirname+'\\'+fname
#
#  To save a path from file open so that you start at that dir next time
#      path=os.path.dirname(FileName)
#      os.chdir(path)
#
#  To write a series of lines into a file:
#    myfile.write(FirstLine+'\n')     note: '\n' is newline char
#    myfile.write(SecondLine+'\n')
#    myfile.write(ThirdLine+'\n')
#
#  Manual file open:  a='C:/MDF-'+sampname+'.txt'
#                     myfile=open(a,'w')
#
#
#  EasyGUI file open:   name=easygui.fileopenbox()  
#                       myfile=open(name,'r')
#  EasyGUI new file:    name=easygui.filesavebox()  
#                       myfile=open(name,'a')
#
#  BETTER: using context manager:  with open(filename,'r') as f:
#                                      thislist=f.readlines()
#
#  Note: context manager automatically
#  closes the file!!!
#
#  NOTE When opening a file...
#  'r' means read - opens file for reading only
#  'a' means append - adds to existing file (or creates new file)
#  'w' means write - makes new file for writing
#                    (vaporizes contents of existing file)
#  
#  If file is open, MyFile.closed is false (and visa versa)
#
#  To close a file:  MyFile.close()
#
#  Error handling
#    try:
#       myfile=open("C:\junk.txt",'r')
#    except:
#       print "No such file!"
#
#  or, if there is a predefined error boolean --
#    try:
#       myfile=open("C:\junk.txt")
#    except IOError:
#       print "No such file!"
#
#  use sys.exit() to blast your way out of a program
#
#  BE CAREFUL passing mutable objects like lists to functions
#  b/c modifications affect object outside scope of function!!
#  Need to remake it in the function and modify the new thing
#
# DISLIN HINTS:
#
#Set the DISLIN environment variable to c:\dislin, the PYTHONPATH
# variable to c:\dislin\python and include c:\dislin\win in your 
# path. 

# The environment variables can be set or modified with the Control Panel 
# (see Control Panel -> System -> Advanced -> Environment Variables).
#
# Note well: Windows fonts dislin.winfnt('ARIAL') looks good on
# screen, but doesn't go to PDF.  You get the default font instead.
# Best option for PDF is dislin.simplx()
#
# This is the order:
#
# 1.  page format, file format, filename
# 2.  dislin.disini()
# 3.  set plot parameters
# 4.  plot axis
# 5.  plot title
# 6.  plot data
# 7.  dislin.disfin()
#
#=============================================================
import math, os, unittest  #from GeoMag.
#apparently no conflict with from import calls that follow

from datetime import date  #from GeoMag.   

from math import sin, radians, cos, atan2, sqrt, degrees, \
                 exp, log, acos, pi, atan, asin, acos, tan

from os import system  #for system("pause") to keep console open at app end

import random   #for fakefish and booter
#==========================================
#Function getdirlist     Mar 2020
#
#Creates a text file list of all files in a dir that you select.

# from easygui import diropenbox

def getdirlist():
    path=diropenbox()
    os.chdir(path)
    file_list=os.listdir(path) 
    myfile=open('DirList.txt','w')
    for i in range(len(file_list)):
        print(file_list[i])
        myfile.write(file_list[i]+'\n')    
    myfile.close()
    print("DirList.txt created.")

#======================================
#Function GetFileNames(SamFileName):   Jan 2016
#
#Give it a SamFileName, and it will return a list of
#file names that you can then open.
#
#If the file you give it does not have .sam extension,
#you will get back a list of one item - that file name.
#
#So, you can also use this to get a single file name for opening also.
#
#Use code like this to get the .sam file name:
#-----------
#  print 'File open window may be hiding beneath other windows.'
#  ans=raw_input('                      Hit enter to open a .sam file: ')
#  FileName=easygui.fileopenbox()
#-------------
#
#Returns a list with the path, samfilename, and list of full file names

def GetFileNames(SamFileName):
    #check to see whether we got a .sam file or data file
    ext=os.path.basename(SamFileName)[-4:]
    if ext in ['.sam','.SAM']:
        samyes=True
    else:
        samyes=False
    if samyes:    
        MyFile=open(SamFileName,'r')
        #Get list of file names.  First 2 will be headers.
        Lines=MyFile.readlines()
        MyFile.close()
        numfiles=len(Lines)
        print()
        print('There are ',numfiles,'lines in the .sam file')
        #get path for sam file using os filename split method
        pathtofile=os.path.dirname(SamFileName)
        samfile=os.path.basename(SamFileName)
        #print 'Folder is '+pathtofile
        #print 'sam file is '+samfile
        #first two lines of sam file are header
        #Use [:-1] to slice off end of line character
        ListOfFiles=list()
        i=2
        filecount=0
        while i<len(Lines):
            #ignore blank lines
            if len(Lines[i])>1:
                fullname=pathtofile+'\\'+Lines[i][:-1]
                ListOfFiles.append(fullname)
                filecount=filecount+1
            i=i+1
        
    else:    # here is we just got a data file rather than a sam file
        print()
        print('You did not pick a sam file.')
        print('The name of the file you opened will be returned.')
        ListOfFiles=[SamFileName]
        pathtofile=os.path.dirname(SamFileName)
        samfile=os.path.basename(SamFileName) 
        filecount=1
            
    print('There are ',filecount,' data files in the returned list')
    print()
    '''
    i=0
    while i<len(ListOfFiles):
        print ListOfFiles[i]      
        i=i+1   
    '''    
    outlist=[pathtofile,samfile,ListOfFiles]
    return outlist
#==============================================
#Function ListFromTextFile()    Oct 2012
#
#Uses easygui to get a file name,
#opens the file and turns the lines
#into a list of text lines.  Closes file.
#Make sure to warn user that file open dialog
#may be hidden!
#
#Modified Jan 2017 to change to directory of first call
#
# from easygui import fileopenbox
def ListFromTextFile():
    FileName=fileopenbox()
    MyFile=open(FileName,'r')
    Lines=MyFile.readlines()
    MyFile.close()
    #get path we have drilled to...
    path=os.path.dirname(FileName)
    #... and start there next time.
    os.chdir(path)
    #print os.getcwd()
    return Lines


#==============================================
#Function ListFromTextFile2()   Mar 2020
#
#Use this one if you want the path to first file
#so that you cwdir to it in calling code
#
#Uses easygui to get a file name,
#opens the file and turns the lines
#into a list of text lines.  Closes file.
#Make sure to warn user that file open dialog
#may be hidden!

#
#Modified Jan 2017 to change to directory of first call
#
# from easygui import fileopenbox
def ListFromTextFile2():
    FileName=fileopenbox()
    MyFile=open(FileName,'r')
    Lines=MyFile.readlines()
    MyFile.close()
    #get path we have drilled to...
    path=os.path.dirname(FileName)
    name=os.path.basename(FileName)
    #... and start there next time.
    # os.chdir(path)
    #print os.getcwd()
    return [Lines,path,name]


#======================================
#Function ListFromTextFileExt(ext)    Oct 2013
#
#Uses easygui to get a file name with a particular extension
#which must be a string, like "rmg" or "txt"
#Opens the file and turns the lines
#into a list of text lines.  Closes file.
#Make sure to warn user that file open dialog
#may be hidden!
#
#Modified Jan 2017 to change to directory of first call
#
# from easygui import fileopenbox
def ListFromTextFileExt(ext):
    FileName=fileopenbox(default='*.'+ext)
    MyFile=open(FileName,'r')
    Lines=MyFile.readlines()
    MyFile.close()
    #get path we have drilled to...
    path=os.path.dirname(FileName)
    #... and start there next time.
    os.chdir(path)
    #print os.getcwd()
    return Lines

#===============================================
#Function ListToTextFile(List of text lines)    Oct 2012
#
#Uses easygui to get a file name,
#opens the file, and then writes (or appends) a list of lines
#to the text file.
#
#Returns the filename   ADDED MAY 5 2014
#
#Make sure that each line in list has
#newline character at the end.  Closes file. 
#
#Make sure to warn user that file open dialog
#may be hidden!
#
# from easygui import filesavebox
def ListToTextFile(ListOfLines):
    FileName=filesavebox()
    MyFile=open(FileName,'a')
    i=0
    while i<len(ListOfLines):
       MyFile.write(ListOfLines[i])
       i=i+1               
    MyFile.close()
    return FileName
    
#===============================================   
#Function getNums
#Gives a list with first N nums in a string.
#Use to extract data from a line in a file.
#If anything goes wrong, puts 99999 in list
#
#
def getNums(textline,n):
  numlist=textline.split()  #splits string into list of words
  outlist=list()          #create an empty list
  i=0
  while i<n:
      try:
         a=float(numlist[i])
      except: #on any error, put 99999 into list
         a=99999
      outlist.append(a)
      i=i+1
  return outlist

#=================
#Function getAllItems  (FOR COMMMA DELIMITED TEXT FILE)
#
#Take a comma delimited text line and returns a list of items in the line.
#Every item that can be turned into a number is converted via float and added to the list.
#Otherwise, item is left as text and added to the list. 
#Use to extract data from a line in a file (comma delimited).
#

def getAllItems(textline):
  numlisttxt=textline.split(',')  #splits string with commas into list of strings
  #How many did we catch?
  n=len(numlisttxt)
  outlist=list()          #create an empty list
  for i in range(n):
      try:
         #does it look like a number?  If so, make it one...
         a=float(numlisttxt[i])
      except:
         #if not, then just put the unconverted text string in the list
         a=numlisttxt[i]
      outlist.append(a)
  return outlist

#=============================================================    
#Function NumListToText
#Converts list of nums to line of spaced text WITH ENDLINE.
#If anything goes wrong, puts '99999' in list
#
def NumListToText(numlist):
  i=0
  textline=''
  while i<len(numlist):
      try:
         a=repr(numlist[i])
      except: #on any error, put 99999 into list
         a='99999'
      textline=textline+a+' '
      i=i+1
  textline=textline+'\n'  #add CR for endline
  return textline

#======================================================
#Function getDIMA
#Gives a list of *TEXT* [Dec,Inc,Mom,error angle,tag] from
#a line from a Oxy/Caltech data file.
#REMEMBER - first two lines in file are header.  Next may be blank.
#
#
def getDIMA(textline):
  outlist=list()          #create an empty list   
  tag=textline[:6]   #first six chars are tag
  data=textline[6:] 
  numlist=data.split()  #splits string into list of words
  outlist.append(numlist[0]) #dec
  outlist.append(numlist[1]) #inc
  outlist.append(numlist[4]) #mom
  outlist.append(numlist[5])  #error angle
  outlist.append(tag)       #tag
  return outlist

#======================================================
#Function getDIMAnums
#Gives a list of 4 numbers plus a string [Dec,Inc,Mom,error angle,tag] from
#a line from a Oxy/Caltech data file.
#REMEMBER - first two lines in file are header.  Next may be blank.
#
#
def getDIMAnums(textline):
  outlist=list()          #create an empty list   
  tag=textline[:6]   #first six chars are tag
  data=textline[6:]  #data will be rest of line of text...
  numlist=getNums(data,6)  #... which should start with 6 nums  
  outlist.append(numlist[0]) #dec
  outlist.append(numlist[1]) #inc
  outlist.append(numlist[4]) #mom
  outlist.append(numlist[5])  #error angle
  outlist.append(tag)       #tag
  return outlist

#======================================================
#Function getstratDIMAnums
#Gives a list of 4 numbers plus a string [stratDec,stratInc,Mom,error angle,tag] from
#a line from a Oxy/Caltech data file.
#REMEMBER - first two lines in file are header.  Next may be blank.
#
#
def getstratDIMAnums(textline):
  outlist=list()          #create an empty list   
  tag=textline[:6]   #first six chars are tag
  data=textline[6:]  #data will be rest of line of text...
  numlist=getNums(data,6)  #... which should start with 6 nums  
  outlist.append(numlist[2]) #dec - tilt corrected
  outlist.append(numlist[3]) #inc - tilt corrected
  outlist.append(numlist[4]) #mom
  outlist.append(numlist[5])  #error angle
  outlist.append(tag)       #tag
  return outlist
#======================================================
#Function getspecDIMAnums
#Gives a list of 4 numbers plus a string [specDec,specInc,Mom,error angle,tag] from
#a line from a Oxy/Caltech data file.
#REMEMBER - first two lines in file are header.  Next may be blank.
#
#
def getspecDIMAnums(textline):
  outlist=list()          #create an empty list   
  tag=textline[:6]   #first six chars are tag
  data=textline[6:]  #data will be rest of line of text...
  numlist=getNums(data,8)  #... which should start with 6 nums  
  outlist.append(numlist[6]) #dec - tilt corrected
  outlist.append(numlist[7]) #inc - tilt corrected
  outlist.append(numlist[4]) #mom
  outlist.append(numlist[5])  #error angle
  outlist.append(tag)       #tag
  return outlist
#======================================================
#Function GetDataFromRMG
#
#Gets treatment,lev,Mz,Dec,Inc,M,ARM bias field, and seconds since last midnight from data line in RMG file
#Input is a list of text lines from the rmg file
#Returns a list of lists.
#First line is samp file name.
#The rest are lists with [treat,lev,Mz,Dec,Inc,M,arm bias field,secs],
#all nums except for treat which is a string
#
#1/2016  Modified GetDataFromRMG to handle old and new RMG formats
#9/2018  Add the time bit to track down funkiness in the IRM system
#
#Lots of new fields added at some point, so old test for data line on number of
#fields no good.
#
#Fixed test for data line so that we avoid those that start with
#something other than a treatment tag.
#
# in RMG: 5=Mz  9=Mx  11=My  14=core dec  15=core inc  16=moment  17=ang  19=date/time
#         0=treatment tag  1=level  2=arm bias field.  Note: old RMGs end at field 12.
#

def GetDataFromRMG(biglist):

    bigoutlist=list()
    #sample name is first item in first line of RMG file
    wordlist=biglist[0].split(',')
    bigoutlist.append(wordlist[0])
    i=1
    while i<len(biglist):
        wordlist=biglist[i].split(',')
        #print wordlist
        #now, look for lines with data, and skip headers
        if wordlist[0]==' ':
            skip=True
        elif wordlist[0][:4] in {'Inst','Time'}:
            skip=True
        else:
            skip=False
            
        if not skip:   #so skip short header line
            outlist=list()
            outlist.append(wordlist[0])  #tag
            
            try:
                a=float(wordlist[1])   #level
            except: #on any error, put 99999 into list
                a=99999
            outlist.append(a)

            try:
                a=float(wordlist[5])   #Mz
            except: #on any error, put 99999 into list
                a=99999
            outlist.append(a)

            try:
                a=float(wordlist[14])   #core dec
            except: #on any error, put 99999 into list
                a=99999
            outlist.append(a)

            try:
                a=float(wordlist[15])   #core inc
            except: #on any error, put 99999 into list
                a=99999
            outlist.append(a)

            try:
                a=float(wordlist[16])  #moment
            except: #on any error, put 99999 into list
                a=99999
            outlist.append(a)

            try:
                a=float(wordlist[2])  #arm bias field
            except: #on any error, put 99999 into list
                a=99999
            outlist.append(a)
            
            try:
                #there might be leading blank...
                if wordlist[19][0] in ['',' ','   ']:
                    wordlist[19]=wordlist[19][1:]
                timeanddate=wordlist[19].split(' ')    #0 is date  1 is time              
                time=timeanddate[1]             
                timelist=time.split(':')
                secs=3600*int(timelist[0])+60*int(timelist[1])+int(timelist[2])
                outlist.append(secs)
            except:
                outlist.append(99999999)

            #if old format RMG, with just 13 fields,
            #then compute core dec, inc, and m from x, y, and z
            if len(wordlist)<20:
                try:
                    dimlist=dim(float(wordlist[9]),float(wordlist[11]),float(wordlist[5]))
                except:
                    #print wordlist[9],wordlist[11],wordlist[5]
                    dimlist=[99999,99999,99999]
                outlist[3]=dimlist[0]   #core dec
                outlist[4]=dimlist[1]   #core inc
                outlist[5]=dimlist[2]   #moment
            
            bigoutlist.append(outlist)
            #print outlist
        i=i+1
        #print i
    return bigoutlist    
#======================================================
#Fisher statistics.
#Input in N and R
#Returns list of [k,alpha-95]
def FishStat(N,R):
    #deal with unlikely input
    if R==0:
        k=0
        a95=180    
    elif R>=N:
        k=99999
        a95=0  
    #almost always here
    else:
        N=float(N)  #make sure its real
        k=((N-1)/(N-R))
        M=pow(20,1/(N-1))
        Q=(N-R)/R
        if (1-(Q*(M-1)))<-1.0: a95=180.0
        else:  a95=degrees(acos(1-(Q*(M-1)))) #usually here
    return[k,a95]

#==============================================================
#Function FakeFish
#Fake sampling of size N of an artificial Fisher distribution
#with mean of dec,inc and precision parameter k.
#Returns a list of N [sampDec,sampInc] pairs.
#Pythonized from the Pascal by Bogue 10-24-05
#
def FakeFish(dec,inc,k,N):
    SampList=list()
    i=0
    while i<N:
        #random() gives real between 0 and 1
        #Oct 2007 changed import to here
        rand1=random.random()
        rand2=random.random()
        lamb=exp(-2*k)
        #note: log with no first arg is ln
        rawCos=sqrt(-log(rand1*(1-lamb)+lamb)/(2*k))
        IRot=2*degrees(acos(rawCos))-90.0  #Fisher I from  90
        DRot=degrees(2*pi*rand2)           #random D
  #Now, rotate center of Fisher dist. to  mean direction
        #and append to list.
        SampList.append(HRot(DRot,IRot,dec-90,90-inc))
        i=i+1
    return SampList

#======================================================
# Fisher function: calculate Fisher mean and stats
# Takes a list of [dec,inc] lists. 
# Returns list of [MeanDec,Meaninc,k,alpha95,R].
# Uses pmag funcs VecSum and FishStat.
# Everything in degrees.
#Hastily coded and lightly tested by Bogue 11-3-05
#Revised 15 Nov 2009 so that passed list not modified
def Fisher(diList):
    #add moment=1 to each d,i pair
    dimList=list()
    i=0
    while i<len(diList):
        thisdim=[]
        thisdi=diList[i]
        thisdim.append(thisdi[0])
        thisdim.append(thisdi[1])
        thisdim.append(1)
        dimList.append(thisdim)
        i=i+1
    DIR=VecSum(dimList)
    kalpha=FishStat(len(dimList),DIR[2])
    outFish=[DIR[0],DIR[1],kalpha[0],kalpha[1],DIR[2]]
    return outFish

#=================================================================
#mixfit: Calculate mean direction from mixed directions and demag circles
#
#McFadden and McElhinny 1988, EPSL 88, 161-172
#
#Takes two lists: dirlist is list of [dec,inc] for directions
#                polelist is list of [dec,inc] for poles to demag circles
#and a boolean:  if just 1 item in dirlist and guess=True, then it's just a guessed
#                 direction to get the iteraion started.  Won't be used after that.
#
#Returns [dec,inc,R,K,a95,N]
#
#test data from McFadden et al 1988. 6 poles and 2 dirs#
# p1=[68.4,-44.4]
# p2=[22.9,-30.5]
# p3=[27.4,-38.8]
# p4=[48.2,-27.2]
# p5=[36.8,-37.3]
# p6=[39.8,-35.2]

# d1=[146.6,-52.6]
# d2=[195.5,-57.9]

#with no "sector constraints":
# 183.1,-51.2  R:7.8336 k:24.032  a95:12.5  N=8
#
#========================================
def mixfit(dirlist,polelist,guess):
    #set direction to project to great circles
    if len(dirlist)>1:
        dirsum=Fisher(dirlist)   #returns D,I,k,a95,R
    else:
        dirsum=[dirlist[0][0],dirlist[0][1]]

    dirsumlist=[dirsum]  #we will append gc points to this in first loop    

    #set up list for direction picked along each great circle. Fill with dummy values for now.    
    gcdirs=list()   
    for i in range(len(polelist)):
        gcdirs.append(polelist[i])

    for i in range (len(polelist)):
        #first run through iteration
        gcdir=closest(dirsum[0],dirsum[1],polelist[i][0],polelist[i][1])
        newdir=[gcdir[0],gcdir[1]]
        gcdirs[i]=newdir  #so replace the dummy value with this first estimate of point on gc
        #update the dirsum direction with the gc dir
        dirsumlist.append(newdir)
        dirsum=Fisher(dirsumlist)
            
    #get rid of first item in dirlist if it is just a seed guess
    if (len(dirlist)==1) and guess:   #need just guess I think
        dirlist=list()    #now empty,so will just use great circles from now on
        
    stop=False
    kk=0
    while not(stop):
        oldgcdirs=list()
        for i in range (len(polelist)):
            #update guess of dir. Start with known dirs, same as before
            newdirsum=list()
            for j in range (len(dirlist)):
                newdirsum.append(dirlist[j])
                
            #delete current element of gcdirs. Do it to a copy of the list
            tempgcdirs=list()
            for j in range(len(gcdirs)):
                tempgcdirs.append(gcdirs[j])
                oldgcdirs.append(gcdirs[j])
            
            del tempgcdirs[i]
            
            #now add remaining gcdirs to the dir guess
            for j in range (len(tempgcdirs)):
                newdirsum.append(tempgcdirs[j])
             
            #and calculate the new guess    
            newguess=Fisher(newdirsum)
            
            newgcdir=closest(newguess[0],newguess[1],polelist[i][0],polelist[i][1])
            #now, replace current element of gcdirs with updated estimate
            gcdirs[i]=[newgcdir[0],newgcdir[1]]

            #lets see how we are doing and stop if estimates are no longer changing much
            angsum=0
            for j in range(len(gcdirs)):            
                angsum=angsum+ang(gcdirs[j][0],gcdirs[j][1],oldgcdirs[j][0],oldgcdirs[j][1])
            if (angsum<0.00001) or (kk>1000):
                stop=True
        kk=kk+1
        #print 'loop #',kk,' angsum is ',angsum
        
    #now get mean direction and stats
    newdirsum=list()
    M=len(dirlist)
    for j in range (M):
        newdirsum.append(dirlist[j])
    N=len(gcdirs)    
    for j in range (N):
        newdirsum.append(gcdirs[j])
    meandir=Fisher(newdirsum)
    #print meandir
    R=meandir[4]
    k=(2.0*M+N-2.0)/(2.0*(M+N-R))
    #print k
    Nprime=M+(N/2.0)
    W=pow(20.0,1.0/(Nprime-1.0))
    Q=(Nprime-1.0)/(k*R)
    if (1-(Q*(W-1.0)))<-1.0:
        a95=180.0
    else:
        a95=degrees(acos(1-(Q*(W-1.0)))) #usually here
    outlist=[meandir[0],meandir[1],R,k,a95,N+M,kk,angsum]
    return outlist

#=========================================================
#function delvgp
#calculates angular standard deviation of VGPs from spin axis
#takes list of [lon,lat]s
#returns asd in degrees   
def delvgp(vgplist):
   i=0
   sumsq=0
   while i<len(vgplist):
      thisvgp=vgplist[i]
      vgpNED=ned(thisvgp[0],thisvgp[1],1.)
      thisang=degrees(acos(vgpNED[2]))  #dot prod with (0,0,1)
      sumsq=sumsq+thisang*thisang
      i=i+1  
   delta=sqrt(sumsq/len(vgplist))
   return(delta)
#=================================================================
#function getasd
#Take a list of lon,lat pairs and gets ASD from its mean and spin axis
#Returns [ASD-mean, ASD-spinaxis ,mean lon, meanl lat]

def getasd(inlist):
    meanpole=Fisher(inlist)
    meanlon=meanpole[0]
    meanlat=meanpole[1]

    summ=0
    sump=0

    N=len(inlist)
    
    for i in range(N):
        a1=ang(inlist[i][0],inlist[i][1],meanlon,meanlat)
        summ=summ+a1*a1
        a2=ang(inlist[i][0],inlist[i][1],90.0,90.0)
        sump=sump+a2*a2
        
    asdm=sqrt(summ/(N-1))
    asdp=sqrt(sump/N)

    return [asdm,asdp,meanlon,meanlat]
              

#======================================================
#DecIncStats function: takes list of [d,i] lists.
#Finds fisher mean, and then rotates distribution
#to horizontal.
#Gives sd of both dec and inc to use as measure of elongation
#that is: elongation=isd/dsd
#Returns [av dec(meand), sd dec(dsd), av inc(meani), sd inc(isd)]
def DecIncStats(dilist):
    #get mean
    Fishout=Fisher(dilist)
    meand=Fishout[0]
    meani=Fishout[1]
    #rotate whole mess to horiz to see dec and inv variance
    rotlist=list()
    i=0
    while i<len(dilist):
        thisdi=dilist[i]
        rotdi=HRot(thisdi[0],thisdi[1],meand-90,meani)
        rotlist.append(rotdi)
        i=i+1
        
#back thru list to get variance of decs
    i=0
    dsum=isum=0
    while i<len(rotlist):
        thisdi=rotlist[i]
        decdiff=thisdi[0]-meand
        if decdiff>180: decdiff=360-decdiff
        elif decdiff<-180: decdiff=360+decdiff
        #sum up squares of diffs
        dsum=dsum+(decdiff*decdiff)
        isum=isum+(thisdi[1])*(thisdi[1])
        i=i+1

    dsd=sqrt(dsum/(len(rotlist)-1))    
    isd=sqrt(isum/(len(rotlist)-1))
    statlist=[meand,dsd,meani,isd]
    return statlist

#====================================================
#
#function booter
#Takes a list and does a bootstrap resampling of it, producing a new list.
#
def booter(inlist):
    bootlist=list()
    n=len(inlist)
    i=0
    while i<n:
        #get random int between 0 and n-1
        bootlist.append(random.choice(inlist))
        i=i+1  
    return (bootlist)

#==============================================================
#PaleoDir function: calculates ancient field from VGP
#Takes (VGPlon,VGPlat,sitelon,sitelat).
#Returns list of [AncDec,AncInc]
#All in degrees.
# Pythonized from Pascal and lightly tested by Bogue 10-28-05
#
def PaleoDir(vlon,vlat,slon,slat):
    if slat>89.95: slat=89.95 #avoids bombouts calculating c below
    sitexyz=ned(slon,slat,1.0)
    vgpxyz=ned(vlon,vlat,1.0)
    sitex=sitexyz[0]
    sitey=sitexyz[1]
    sitez=sitexyz[2]
    vgpx=vgpxyz[0]
    vgpy=vgpxyz[1]
    vgpz=vgpxyz[2]
    delta=acos(sitex*vgpx+sitey*vgpy+sitez*vgpz) #angle pole to site
    a=           cos(radians(90.-vlat))
    b=cos(delta)*cos(radians(90.-slat))
    c=sin(delta)*sin(radians(90.-slat))
    if c>(a-b): dec=degrees(acos((a-b)/c))
    else: dec=0.0   #watch out for bombouts from rounding errors
    if (vlon-slon)>180 and (vlon-slon)<360: dec=360.-dec
    if (vlon-slon)>-180 and (vlon-slon)<0: dec=360.-dec
    if delta==0:inc=90  #another bomb site
    else: inc=degrees(atan(2*cos(delta)/sin(delta)))
    return[dec,inc]

#======================================================
# VGP function: calculate VGP position
# Takes dec,inc at site defined by lon,lat
# Returns list of [VGPlon,VGPlat].
# Everything in degrees.
# Pythonized from Pascal and lightly tested by Bogue 10-25-05
def vgp(dec,inc,lon,lat):
    if inc==0: inc=0.0000000001
    incR=radians(inc)
    latR=radians(lat)
    decR=radians(dec)
    P=atan(2*cos(incR)/sin(incR))
    if P<0: P=P+pi
    VGPlat=asin(sin(latR)*cos(P)+cos(latR)*sin(P)*cos(decR))
    if VGPlat>=pi/2: VGPlat=VGPlat-0.00000001
    if VGPlat<=-pi/2:VGPlat=VGPlat+0.00000001
    B=asin(sin(P)*sin(decR)/cos(VGPlat))
    if cos(P)>=sin(latR)*sin(VGPlat):
         VGPlon=lon+degrees(B)
    else:
         VGPlon=lon+180-degrees(B)
    if VGPlon<0: VGPlon=VGPlon+360
    if VGPlon>=360: VGPlon=VGPlon-360
    vgpList=[VGPlon,degrees(VGPlat)]    
    return vgpList
#======================================================
#vgperr function: calculate VGP with error ellipse 
#from dec,inc, site coords, and and alpha-95 on direction
#dp=semiaxis along site-vgp great circle  dm=other semiaxis
#b95=radius of error circle with same area as dp-dm ellipse
#Everything in degrees
#Returns list [VGPlon,VGPlat,dp,dm,b95]
#lightly tested aug 2013
def vgperr(dec,inc,slon,slat,a95):
    plist=vgp(dec,inc,slon,slat)
    #print plist
    p=atan2(2,tan(radians(inc)))
    #print degrees(p)  
    dp=2*a95*(1/(1+3*cos(radians(inc))*cos(radians(inc))))
    #print dp
    dm=a95*(sin(p)/cos(radians(inc)))
    #print dm
    b95=sqrt(dm*dp)
    plist.append(dp)
    plist.append(dm)
    plist.append(b95)
    return plist

#=====================================================
#modelG-k function: calculate k for VGPs according to Model G
#For 0-5 Ma
def kmodG(lat):
    ssq=(0.23*lat)*(0.23*lat)+12.8*12.8
    k=6561/ssq    #n.b.  s=sqrt(2/k) with s in radians
    return k
#=====================================================
#modelG-k function: calculate k for VGPs according to Model G
#For Miocene 5-22.5 Ma
def kmodGmio(lat):
    ssq=(0.19*lat)*(0.19*lat)+17.8*17.8
    k=6561/ssq
    return k
#========================================================
#Flatten I function.
#Takes I(H) and f, give I(DRM).  All degrees
def Flatten(Ih,f):
    Ih=radians(Ih)
    flatinc=degrees(atan(f*tan(Ih)))
    return(flatinc)

#========================================================
#UnFlatten I function.
#Takes I(DRM) and f, give I(H).  All degrees
def UnFlatten(Idrm,f):
    Idrm=radians(Idrm)
    Unflatinc=degrees(atan((1/f)*tan(Idrm)))
    return(Unflatinc)

#===========================================================
#GeoMag Class
#calculates field from WMM10
#Instantiate it like this:
#
#WMM10=GeoMag("WMM.COF")
#NDF10=GeoMag("WMM-NDF.COF")
#NAD10=GeoMag("WMM-NAD.COF")
#
#Then, calls like this for field at some lat,lon, elev=0, on 1/1/2010
#
#total = WMM10.GeoMag(Lat, Lon,0,date(2010,1,1))
#NDF = NDF10.GeoMag(Lat,Lon,0,date(2010,1,1))
#NAD = NAD10.GeoMag(Lat,Lon,0,date(2010,1,1))
#
#The answers are: total.dec, total.dip, total.ti, total.bx, total.by, total.bz
#
#Needs a text file of cofficients to call
#WMM.COF is the basic one
#WMM-NAD.COF and WMM-NDF.COF are with axial dipole and all dipole terms set to 0
#
#WMM.COF looks like this:
'''
    2010.0            WMM-2010        11/20/2009
  1  0  -29496.6       0.0       11.6        0.0
  1  1   -1586.3    4944.4       16.5      -25.9
  2  0   -2396.6       0.0      -12.1        0.0
  2  1    3026.1   -2707.7       -4.4      -22.5
  2  2    1668.6    -576.1        1.9      -11.8
  3  0    1340.1       0.0        0.4        0.0
  3  1   -2326.2    -160.2       -4.1        7.3
  3  2    1231.9     251.9       -2.9       -3.9
  3  3     634.0    -536.6       -7.7       -2.6
  4  0     912.6       0.0       -1.8        0.0
  4  1     808.9     286.4        2.3        1.1
  4  2     166.7    -211.2       -8.7        2.7
  4  3    -357.1     164.3        4.6        3.9
  4  4      89.4    -309.1       -2.1       -0.8
  5  0    -230.9       0.0       -1.0        0.0
  5  1     357.2      44.6        0.6        0.4
  5  2     200.3     188.9       -1.8        1.8
  5  3    -141.1    -118.2       -1.0        1.2
  5  4    -163.0       0.0        0.9        4.0
  5  5      -7.8     100.9        1.0       -0.6
  6  0      72.8       0.0       -0.2        0.0
  6  1      68.6     -20.8       -0.2       -0.2
  6  2      76.0      44.1       -0.1       -2.1
  6  3    -141.4      61.5        2.0       -0.4
  6  4     -22.8     -66.3       -1.7       -0.6
  6  5      13.2       3.1       -0.3        0.5
  6  6     -77.9      55.0        1.7        0.9
  7  0      80.5       0.0        0.1        0.0
  7  1     -75.1     -57.9       -0.1        0.7
  7  2      -4.7     -21.1       -0.6        0.3
  7  3      45.3       6.5        1.3       -0.1
  7  4      13.9      24.9        0.4       -0.1
  7  5      10.4       7.0        0.3       -0.8
  7  6       1.7     -27.7       -0.7       -0.3
  7  7       4.9      -3.3        0.6        0.3
  8  0      24.4       0.0       -0.1        0.0
  8  1       8.1      11.0        0.1       -0.1
  8  2     -14.5     -20.0       -0.6        0.2
  8  3      -5.6      11.9        0.2        0.4
  8  4     -19.3     -17.4       -0.2        0.4
  8  5      11.5      16.7        0.3        0.1
  8  6      10.9       7.0        0.3       -0.1
  8  7     -14.1     -10.8       -0.6        0.4
  8  8      -3.7       1.7        0.2        0.3
  9  0       5.4       0.0       -0.0        0.0
  9  1       9.4     -20.5       -0.1       -0.0
  9  2       3.4      11.5        0.0       -0.2
  9  3      -5.2      12.8        0.3        0.0
  9  4       3.1      -7.2       -0.4       -0.1
  9  5     -12.4      -7.4       -0.3        0.1
  9  6      -0.7       8.0        0.1       -0.0
  9  7       8.4       2.1       -0.1       -0.2
  9  8      -8.5      -6.1       -0.4        0.3
  9  9     -10.1       7.0       -0.2        0.2
 10  0      -2.0       0.0        0.0        0.0
 10  1      -6.3       2.8       -0.0        0.1
 10  2       0.9      -0.1       -0.1       -0.1
 10  3      -1.1       4.7        0.2        0.0
 10  4      -0.2       4.4       -0.0       -0.1
 10  5       2.5      -7.2       -0.1       -0.1
 10  6      -0.3      -1.0       -0.2       -0.0
 10  7       2.2      -3.9        0.0       -0.1
 10  8       3.1      -2.0       -0.1       -0.2
 10  9      -1.0      -2.0       -0.2        0.0
 10 10      -2.8      -8.3       -0.2       -0.1
 11  0       3.0       0.0        0.0        0.0
 11  1      -1.5       0.2        0.0       -0.0
 11  2      -2.1       1.7       -0.0        0.1
 11  3       1.7      -0.6        0.1        0.0
 11  4      -0.5      -1.8       -0.0        0.1
 11  5       0.5       0.9        0.0        0.0
 11  6      -0.8      -0.4       -0.0        0.1
 11  7       0.4      -2.5       -0.0        0.0
 11  8       1.8      -1.3       -0.0       -0.1
 11  9       0.1      -2.1        0.0       -0.1
 11 10       0.7      -1.9       -0.1       -0.0
 11 11       3.8      -1.8       -0.0       -0.1
 12  0      -2.2       0.0       -0.0        0.0
 12  1      -0.2      -0.9        0.0       -0.0
 12  2       0.3       0.3        0.1        0.0
 12  3       1.0       2.1        0.1       -0.0
 12  4      -0.6      -2.5       -0.1        0.0
 12  5       0.9       0.5       -0.0       -0.0
 12  6      -0.1       0.6        0.0        0.1
 12  7       0.5      -0.0        0.0        0.0
 12  8      -0.4       0.1       -0.0        0.0
 12  9      -0.4       0.3        0.0       -0.0
 12 10       0.2      -0.9        0.0       -0.0
 12 11      -0.8      -0.2       -0.1        0.0
 12 12       0.0       0.9        0.1        0.0
999999999999999999999999999999999999999999999999
999999999999999999999999999999999999999999999999
'''
class GeoMag:

    def GeoMag(self, dlat, dlon, h=0, time=date.today()): # latitude (decimal degrees), longitude (decimal degrees), altitude (feet), date
        #time = date('Y') + date('z')/365
        time = time.year+((time - date(time.year,1,1)).days/365.0)
        alt = h/3280.8399

        otime = oalt = olat = olon = -1000.0

        dt = time - self.epoch
        glat = dlat
        glon = dlon
        rlat = math.radians(glat)
        rlon = math.radians(glon)
        srlon = math.sin(rlon)
        srlat = math.sin(rlat)
        crlon = math.cos(rlon)
        crlat = math.cos(rlat)
        srlat2 = srlat*srlat
        crlat2 = crlat*crlat
        self.sp[1] = srlon
        self.cp[1] = crlon

        #/* CONVERT FROM GEODETIC COORDS. TO SPHERICAL COORDS. */
        if (alt != oalt or glat != olat):
            q = math.sqrt(self.a2-self.c2*srlat2)
            q1 = alt*q
            q2 = ((q1+self.a2)/(q1+self.b2))*((q1+self.a2)/(q1+self.b2))
            ct = srlat/math.sqrt(q2*crlat2+srlat2)
            st = math.sqrt(1.0-(ct*ct))
            r2 = (alt*alt)+2.0*q1+(self.a4-self.c4*srlat2)/(q*q)
            r = math.sqrt(r2)
            d = math.sqrt(self.a2*crlat2+self.b2*srlat2)
            ca = (alt+d)/r
            sa = self.c2*crlat*srlat/(r*d)

        if (glon != olon):
            for m in range(2,self.maxord+1):
                self.sp[m] = self.sp[1]*self.cp[m-1]+self.cp[1]*self.sp[m-1]
                self.cp[m] = self.cp[1]*self.cp[m-1]-self.sp[1]*self.sp[m-1]

        aor = self.re/r
        ar = aor*aor
        br = bt = bp = bpp = 0.0
        for n in range(1,self.maxord+1):
            ar = ar*aor
            
            #for (m=0,D3=1,D4=(n+m+D3)/D3;D4>0;D4--,m+=D3):
            m=0
            D3=1
            #D4=(n+m+D3)/D3
            D4=(n+m+1)
            while D4>0:

        # /*
                # COMPUTE UNNORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
                # AND DERIVATIVES VIA RECURSION RELATIONS
        # */
                if (alt != oalt or glat != olat):
                    if (n == m):
                        self.p[m][n] = st * self.p[m-1][n-1]
                        self.dp[m][n] = st*self.dp[m-1][n-1]+ct*self.p[m-1][n-1]

                    elif (n == 1 and m == 0):
                        self.p[m][n] = ct*self.p[m][n-1]
                        self.dp[m][n] = ct*self.dp[m][n-1]-st*self.p[m][n-1]

                    elif (n > 1 and n != m):
                        if (m > n-2):
                            self.p[m][n-2] = 0
                        if (m > n-2):
                            self.dp[m][n-2] = 0.0
                        self.p[m][n] = ct*self.p[m][n-1]-self.k[m][n]*self.p[m][n-2]
                        self.dp[m][n] = ct*self.dp[m][n-1] - st*self.p[m][n-1]-self.k[m][n]*self.dp[m][n-2]

        # /*
                # TIME ADJUST THE GAUSS COEFFICIENTS
        # */
                if (time != otime):
                    self.tc[m][n] = self.c[m][n]+dt*self.cd[m][n]
                    if (m != 0):
                        self.tc[n][m-1] = self.c[n][m-1]+dt*self.cd[n][m-1]

        # /*
                # ACCUMULATE TERMS OF THE SPHERICAL HARMONIC EXPANSIONS
        # */
                par = ar*self.p[m][n]
                
                if (m == 0):
                    temp1 = self.tc[m][n]*self.cp[m]
                    temp2 = self.tc[m][n]*self.sp[m]
                else:
                    temp1 = self.tc[m][n]*self.cp[m]+self.tc[n][m-1]*self.sp[m]
                    temp2 = self.tc[m][n]*self.sp[m]-self.tc[n][m-1]*self.cp[m]

                bt = bt-ar*temp1*self.dp[m][n]
                bp = bp + (self.fm[m] * temp2 * par)
                br = br + (self.fn[n] * temp1 * par)
        # /*
                    # SPECIAL CASE:  NORTH/SOUTH GEOGRAPHIC POLES
        # */
                if (st == 0.0 and m == 1):
                    if (n == 1):
                        self.pp[n] = self.pp[n-1]
                    else:
                        self.pp[n] = ct*self.pp[n-1]-self.k[m][n]*self.pp[n-2]
                    parp = ar*self.pp[n]
                    bpp = bpp + (self.fm[m]*temp2*parp)
                    
                D4=D4-1
                m=m+1

        if (st == 0.0):
            bp = bpp
        else:
            bp = bp/st
        # /*
            # ROTATE MAGNETIC VECTOR COMPONENTS FROM SPHERICAL TO
            # GEODETIC COORDINATES
        # */
        bx = -bt*ca-br*sa
        by = bp
        bz = bt*sa-br*ca
        # /*
            # COMPUTE DECLINATION (DEC), INCLINATION (DIP) AND
            # TOTAL INTENSITY (TI)
        # */
        bh = math.sqrt((bx*bx)+(by*by))
        ti = math.sqrt((bh*bh)+(bz*bz))
        dec = math.degrees(math.atan2(by,bx))
        dip = math.degrees(math.atan2(bz,bh))
        # /*
            # COMPUTE MAGNETIC GRID VARIATION IF THE CURRENT
            # GEODETIC POSITION IS IN THE ARCTIC OR ANTARCTIC
            # (I.E. GLAT > +55 DEGREES OR GLAT < -55 DEGREES)

            # OTHERWISE, SET MAGNETIC GRID VARIATION TO -999.0
        # */
        gv = -999.0
        if (math.fabs(glat) >= 55.):
            if (glat > 0.0 and glon >= 0.0):
                gv = dec-glon
            if (glat > 0.0 and glon < 0.0):
                gv = dec+math.fabs(glon);
            if (glat < 0.0 and glon >= 0.0):
                gv = dec+glon
            if (glat < 0.0 and glon < 0.0):
                gv = dec-math.fabs(glon)
            if (gv > +180.0):
                gv = gv - 360.0
            if (gv < -180.0):
                gv = gv + 360.0

        otime = time
        oalt = alt
        olat = glat
        olon = glon

        class RetObj:
            pass
        retobj = RetObj()
        retobj.dec = dec
        retobj.dip = dip
        retobj.ti = ti
        retobj.bh = bh
        retobj.bx = bx
        retobj.by = by
        retobj.bz = bz
        retobj.lat = dlat
        retobj.lon = dlon
        retobj.alt = h
        retobj.time = time

        return retobj

    def __init__(self, wmm_filename=None):
        if not wmm_filename:
            wmm_filename = os.path.join(os.path.dirname(__file__), 'WMM.COF')
        wmm=[]
        with open(wmm_filename) as wmm_file:
            for line in wmm_file:
                linevals = line.strip().split()
                if len(linevals) == 3:
                    self.epoch = float(linevals[0])
                    self.model = linevals[1]
                    self.modeldate = linevals[2]
                elif len(linevals) == 6:
                    linedict = {'n': int(float(linevals[0])),
                    'm': int(float(linevals[1])),
                    'gnm': float(linevals[2]),
                    'hnm': float(linevals[3]),
                    'dgnm': float(linevals[4]),
                    'dhnm': float(linevals[5])}
                    wmm.append(linedict)

        z = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        self.maxord = self.maxdeg = 12
        self.tc = [z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13]]
        self.sp = z[0:14]
        self.cp = z[0:14]
        self.cp[0] = 1.0
        self.pp = z[0:13]
        self.pp[0] = 1.0
        self.p = [z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14]]
        self.p[0][0] = 1.0
        self.dp = [z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13]]
        self.a = 6378.137
        self.b = 6356.7523142
        self.re = 6371.2
        self.a2 = self.a*self.a
        self.b2 = self.b*self.b
        self.c2 = self.a2-self.b2
        self.a4 = self.a2*self.a2
        self.b4 = self.b2*self.b2
        self.c4 = self.a4 - self.b4

        self.c = [z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14]]
        self.cd = [z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14]]
        
        for wmmnm in wmm:
            m = wmmnm['m']
            n = wmmnm['n']
            gnm = wmmnm['gnm']
            hnm = wmmnm['hnm']
            dgnm = wmmnm['dgnm']
            dhnm = wmmnm['dhnm']
            if (m <= n):
                self.c[m][n] = gnm
                self.cd[m][n] = dgnm
                if (m != 0):
                    self.c[n][m-1] = hnm
                    self.cd[n][m-1] = dhnm

        #/* CONVERT SCHMIDT NORMALIZED GAUSS COEFFICIENTS TO UNNORMALIZED */
        self.snorm = [z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13]]
        self.snorm[0][0] = 1.0
        self.k = [z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13]]
        self.k[1][1] = 0.0
        self.fn = [0.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0]
        self.fm = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0]
        for n in range(1,self.maxord+1):
            self.snorm[0][n] = self.snorm[0][n-1]*(2.0*n-1)/n
            j=2.0
            #for (m=0,D1=1,D2=(n-m+D1)/D1;D2>0;D2--,m+=D1):
            m=0
            D1=1
            D2=(n-m+D1)/D1
            while (D2 > 0):
                self.k[m][n] = (((n-1)*(n-1))-(m*m))/((2.0*n-1)*(2.0*n-3.0))
                if (m > 0):
                    flnmj = ((n-m+1.0)*j)/(n+m)
                    self.snorm[m][n] = self.snorm[m-1][n]*math.sqrt(flnmj)
                    j = 1.0
                    self.c[n][m-1] = self.snorm[m][n]*self.c[n][m-1]
                    self.cd[n][m-1] = self.snorm[m][n]*self.cd[n][m-1]
                self.c[m][n] = self.snorm[m][n]*self.c[m][n]
                self.cd[m][n] = self.snorm[m][n]*self.cd[m][n]
                D2=D2-1
                m=m+D1

class GeoMagTest(unittest.TestCase):

    d1=date(2010,1,1)
    d2=date(2012,7,1)
    
    test_values = (
        # date, alt, lat, lon, var
        (d1, 0, 80, 0, -6.13),
        (d1, 0, 0, 120, 0.97),
        (d1, 0, -80, 240, 70.21),
        (d1, 328083.99, 80, 0, -6.57),
        (d1, 328083.99, 0, 120, 0.94),
        (d1, 328083.99, -80, 240, 69.62),
        (d2, 0, 80, 0, -5.21),
        (d2, 0, 0, 120, 0.88),
        (d2, 0, -80, 240, 70.04),
        (d2, 328083.99, 80, 0, -5.63),
        (d2, 328083.99, 0, 120, 0.86),
        (d2, 328083.99, -80, 240, 69.45),
    )
    
    def test_declination(self):
        gm = GeoMag()
        for values in self.test_values:
            calcval=gm.GeoMag(values[2], values[3], values[1], values[0])
            self.assertAlmostEqual(values[4], calcval.dec, 2, 'Expected %s, result %s' % (values[4], calcval.dec))

#if __name__ == '__main__':
#    unittest.main()

#===================================================================
#GetField
#Uses GeoMag to get field from WMM-style cofficients
#Needs a file of coefflcients like WMM.COF
#pass file as model.

#Takes (lon,lat,model)
#Returns (Dec,Inc,H,x,y,z)

def GetField(Lon,Lat,model):    
    WMM=GeoMag(model)
    field = WMM.GeoMag(Lat, Lon,0,date(2010,1,1))       
    dec=field.dec  #get rid of minus declinations
    if dec<0:
        dec=dec+360.0
    outlist=(dec,field.dip,field.ti,field.bx,field.by,field.bz)
    return outlist
#===================================================================

#GetField2010
#Uses GeoMag to get field from WMM cofficients
#Needs a file of coefflcients like WMM.COF
#model paramemter:
#    1=WMM10 "WMM.COF"
#    2=NAD  WMM with axial dipole=0  "WMM-NAD.COF"
#    3=NDF  WMM with all dipole=0    "WMM-NDF.COF"

#Takes (lon,lat,model)
#Returns (Dec,Inc,H,x,y,z)
def GetField2010(Lon,Lat,model):
    
    if model==1:
        WMM10=GeoMag("WMM2010-real-one.COF")
        field = WMM10.GeoMag(Lat, Lon,0,date(2010,1,1))
        
    if model==2:
        NAD10=GeoMag("WMM-NAD.COF")
        field = NAD10.GeoMag(Lat, Lon,0,date(2010,1,1))  

    if model==3:      
        NDF10=GeoMag("WMM-NDF.COF")
        field = NDF10.GeoMag(Lat, Lon,0,date(2010,1,1))
        
    dec=field.dec  #get rid of minus declinations
    if dec<0:
        dec=dec+360.0
    outlist=(dec,field.dip,field.ti,field.bx,field.by,field.bz)
    return outlist
#===================================================================


#HRot  - Rotation about horizontal axis.
#CW Rotation of rot degrees about horizontal axis
#whose end is specified by azm.
#Takes (dec,inc,azm, rot) all in degrees.
#Returns list of [NewDec,NewInc] in degrees
#Pythonized from Pascal and lightly tested by Bogue 10-24-05
#
#To do a strike/dip correction, use strike for azm and rot for the dip.
#
def HRot(dec,inc,azm,rot):
    #calc some trigs
    SinA=sin(radians(azm))
    SinR=sin(radians(rot))
    CosA=cos(radians(azm))
    CosR=cos(radians(rot))
    #polar to rec
    nedList=ned(dec,inc,1)  #seems we can use other funcs in pmag!
    N=nedList[0]
    E=nedList[1]
    D=nedList[2]
    #Do the rotation
    #Note:  backslash continues statement on next line
    NewN= N*(SinA*SinA*CosR + CosA*CosA)     \
         +E*(SinA*CosA * (1.0 - CosR))       \
         -D*(SinA*SinR)
    NewE= N*(SinA*CosA * (1.0 - CosR))       \
         +E*(SinA*SinA + CosA*CosA*CosR)     \
         +D*(CosA*SinR)
    NewD= N*(SinA*SinR)                      \
         -E*(CosA*SinR)                      \
         + D*CosR;
    #rec to polar
    NewDIM=dim(NewN,NewE,NewD)
    #finally, make list to return
    NewDec=NewDIM[0]
    NewInc=NewDIM[1]
    NewDI=[NewDec,NewInc]
    return NewDI

#===============================================
#ARot: Rotation about arbitrary axis.
#Brutish CW rotation of rot degrees about axis whose
#end is specified by azm and plg (azimuth and plunge) in degrees.
#Input is (dec,inc,azm,plg,rot)
#Returns list of [NewDec,NewInc] in degrees.
#Not pretty; not fast; usually runs.
#Pythonized from Pascal and lightly tested by Bogue 10-24-05
def ARot(dec,inc,azm,plg,rot):
    #first, rotate the axis to horizontal, bring D,I along 
    azmComp=azm-90.0
    NewDI=HRot(dec,inc,azmComp,plg)
    #now, do rotation of D,I about horiz. axis
    NewDI=HRot(NewDI[0],NewDI[1],azm,rot)
    #finally, rotate axis back, bring D,I along
    NewDI=HRot(NewDI[0],NewDI[1],azmComp,-plg)
    return NewDI

#=====================================================
#doyz function: do USGS (and CIT) rotation from spec to
#geographic coords.   Y and Z are strike of core top.
def doyz(dec,inc,y,z):
    newdi=HRot(dec,inc,270.0,z)
    newdec=newdi[0]+y-90.0
    if newdec<0.0:
        newdec=newdec+360.0
    if newdec>360.0:
        newdec=newdec-360.0
    newdi[0]=newdec
    return newdi

#======================================================   
#DoYZ   Applies Pomeroy YZ corr to spec coord dir
def DoYZ(d,i,y,z):
   newdi=HRot(d,i,270,z)
   gdec=newdi[0]+y-90
   ginc=newdi[1]
   if gdec>360:
      gdec=gdec-360.0
   if gdec<0:
      gdec=gdec+360.0
   return(gdec,ginc)
#======================================================   
#undoyz  Un-Applies Pomeroy YZ corr to geographic coord dir
def undoyz(d,i,y,z):
   cdec=d+90.0-y
   newdi=HRot(cdec,i,270.0,-z)
   cinc=newdi[1]
   cdec=newdi[0]
   if cdec>360:
      cdec=cdec-360.0
   if cdec<0:
      cdec=cdec+360.0
   return(cdec,cinc)
#=================================
#Spherical to rectangular function.
#Takes Dec and Inc (in degrees) and Mom.
#Returns list of [north,east,down]
#==================================
def ned(D,I,M):
   down=M*sin(radians(I))
   Horiz=M*cos(radians(I))
   north=Horiz*cos(radians(D))
   east=Horiz*sin(radians(D))
   nedlist=[north,east,down]
   return nedlist

#====================================
#Rectangular to spherical function
#Takes north, east, and down comps.
#Returns list of [dec,inc,mom]
#=================================
def dim (N,E,D):
   dec=degrees(atan2(E,N))
   if dec<360:  dec=dec+360
   if dec>360: dec=dec-360  #?? seem to need this
   Horiz=sqrt(N*N+E*E)
   inc=degrees(atan2(D,Horiz))
   mom=sqrt(N*N+E*E+D*D)
   dimlist=[dec,inc,mom]
   return dimlist

#=====================================
#Rect to spherical, but takes [D,I,M] as list
def nedlist(dimlist):
    outlist=ned(dimlist[0],dimlist[1],dimlist[2])
    return outlist
#====================================
#Spherical to rect, but takes [N,E,D] as list
def dimlist(nedlist):
    outlist=dim(nedlist[0],nedlist[1],nedlist[2])
    return outlist
#====================================
#veclen.  Takes ned list and gives length
def veclen(vec1):
   vlen=sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2])
   return vlen

#====================================
#Dot product.  Takes  two ned lists and
#gives list with a dot b, ang a-b
#=================================
def dot (vec1,vec2):
   dotnum=vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2]
   arg=dotnum/(veclen(vec1)*veclen(vec2))
   #look out for errors from values just on the far side of 1 or -1
   if arg>1.0:
      arg=1.0
   if arg<-1.0:
      arg=-1.0
   angle=acos(arg)
   angle=angle*180/pi  #degrees please
   dotlist=[dotnum,angle]
   return dotlist
#===================================
#angle between two unit vectors
#just give it decs and incs; returns the angle in degrees
def ang(d1,i1,d2,i2):
    ned1=ned(d1,i1,1.0)
    ned2=ned(d2,i2,1.0)
    dotlist=dot(ned1,ned2)
    return dotlist[1]
#====================================
#Cross product.  Takes  two ned lists and
#gives list with ned of a x b
#=================================
def cross (vec1,vec2):
   n=vec1[1]*vec2[2]-vec1[2]*vec2[1]
   e=vec1[2]*vec2[0]-vec1[0]*vec2[2]
   d=vec1[0]*vec2[1]-vec1[1]*vec2[0]
   crosslist=[n,e,d]
   return crosslist

#====================================
#Cross product as unit vector.  Takes  two ned lists and
#gives list with ned of a x b with length 1
#=================================
def crossunit (vec1,vec2):
   ucrosslist=cross(vec1,vec2)
   xlen=veclen(ucrosslist)
   ucrosslist=[ucrosslist[0]/xlen,ucrosslist[1]/xlen,ucrosslist[2]/xlen]
   return ucrosslist

#===============================================
#Vector Sum function
#Takes a list of vectors and gives their vector sum.
#Each item in list is a list of (dec,inc,mom)
#Returns list of [MeanDec,MeanInc,R]
def VecSum(VecList):
    numlines=len(VecList)
    Nsum=Esum=Zsum=0
    i=0
    while i<numlines:  
        Numlist=VecList[i]
        nedlist=ned(Numlist[0],Numlist[1],Numlist[2])
        Nsum=Nsum+nedlist[0]
        Esum=Esum+nedlist[1]
        Zsum=Zsum+nedlist[2]
        i=i+1
    dimlist=dim(Nsum, Esum, Zsum)
    dec=dimlist[0]
    inc=dimlist[1]
    R=dimlist[2]
    if dec<0: dec=dec+360
    return[dec,inc,R]
   
#==========================================================   
#Decompose a vector along two other directions.
#Gives moms along d1 and d2, and a ratio that gets big if vector is out of plane
#Input vector is d,i,m   dir 1 is d1,i1   dir 2 is d2,i2
def decomp(d,i,m,d1,i1,d2,i2):

    degtorad=pi/180

    #input vector
    inned=ned(d,i,m)
    vec1ned=ned(d1,i1,1)
    vec2ned=ned(d2,i2,1) 

    #find unit vector perpendicular to d1 and d2
    culist=crossunit(vec1ned,vec2ned)  #ned
    poledim=dim(culist[0],culist[1],culist[2])
    #print culist
    #print dim(culist[0],culist[1],culist[2])
    #print poledim

    #now, project input vector onto this unit vector to see how much of it
    #is out of the the plane of components
    #Note: inproj should be small compared to m if this decomp makes sense
    inproj=dot(inned,culist)
    #print(inproj)  #inproj[1] is angle from vec to normal

    #Now, find out how long is the input vector projected on the comp plane
    momproj=sqrt(m*m-inproj[0]*inproj[0])
    #print momproj

    #Find vector perp to normal and input vector.  In x norm gives vector CW of in
    RotAxis=crossunit(inned,culist)
    #print RotAxis
    RotAxisDIM=dim(RotAxis[0],RotAxis[1],RotAxis[2])

    #now, rotate input vector down to the comp plane using 90-inproj[1]
    #note: Use momproj for length of this projected input vector
    InVecPl=ARot(d,i,RotAxisDIM[0],RotAxisDIM[1],90-inproj[1])

    #Now, rotate all three coplanar vectors to horizontal for checking
    RotAzm=poledim[0]+90
    Rot=90-poledim[1]
    invecH=HRot(InVecPl[0],InVecPl[1],RotAzm,Rot) #di
    #print invecH
    dir1H=HRot(d1,i1,RotAzm,Rot)  #di
    dir2H=HRot(d2,i2,RotAzm,Rot)  #di
    #print dir1H
    #print dir2H
    invecHned=ned(invecH[0],invecH[1],1)
    dir1Hned=ned(dir1H[0],dir1H[1],1)
    dir2Hned=ned(dir2H[0],dir2H[1],1)
    check1=crossunit(dir1Hned,invecHned)
    #print check1
    check2=crossunit(invecHned,dir2Hned)
    #print check2

    DecompGood=False
    nodir1=False
    nodir2=False
    nodir12=False
    
    if check1[2]>0 and check2[2]>0:
       DecompGood=True
   
    if DecompGood:
       print("Decomp makes sense")
    else:
       print("Decomp is no good")

    if check1[2]<0 and check2[2]>0:
       nodir2=True  #input has neg dir2 in it

    if check1[2]>0 and check2[2]<0:
       nodir1=True  #input has neg dir1 in it

    if check1[2]<0 and check2[2]<0:
       nodir12=True  #inpupt has neg dir1 and dir2 in it
  
    #Now find angle between proj invec and comp1 and comp2
    InVecPlNED=ned(InVecPl[0],InVecPl[1],1)
    DotList=dot(InVecPlNED,vec1ned)
    Ang1=DotList[1]
    DotList=dot(InVecPlNED,vec2ned)
    Ang2=DotList[1]
    #print Ang1,Ang2
    Ang3=180-(Ang1+Ang2)

    #Now, use law of sines to find the components
    #SEEMS TO BE A PROBLEM WITH SIGN OF MOM - check if comp is
    #in dir of d1(d2)?  sometimes negative?

    if (Ang3==180) or (Ang3==0):
        Mom2=momproj  #minus sign sometimes? Is this right?
    else:
        Mom2=sin(Ang1*degtorad)*momproj/sin(Ang3*degtorad)
    
    if (Ang1==0) or (Ang1==180):
        Mom1=momproj  #minus sign sometimes?
    else:
        Mom1=sin(Ang2*degtorad)*Mom2/sin(Ang1*degtorad)

    projQ=inproj[0]/momproj
    #print Mom1,Mom2,projQ
    if nodir1:
       Mom1=0
    if nodir2:
       Mom2=0
    if nodir12:
       Mom1=0
       Mom2=0
    decomplist=[Mom1,Mom2,projQ,DecompGood]
    return decomplist

#============================================
#localtogeo: Use this to take a horizontal vector at earth's
#surface and recast it as a vector from earth's center
#
#Wants lon and lat of site on surface, dec and inc of direction there.
#
def localtogeo(sitelon,sitelat,dec,inc):
    #set up the rotation matrix
    #Cols are (local) NED axes in geocentric coords
    lam=radians(sitelon)   #spot on surface where we want local vector
    phi=radians(sitelat)
    r1=-sin(phi)*cos(lam)
    r2=-sin(lam)
    r3=-cos(phi)*cos(lam)
    r4=-sin(phi)*sin(lam)
    r5=cos(lam)
    r6=-cos(phi)*sin(lam)
    r7=cos(phi)
    r8=0.0
    r9=-sin(phi)
    
    #pole  to rec the site where we want local direction
    sitevec=ned(sitelon,sitelat,1.0)

    #pole to rec the direction
    localvec=ned(dec,inc,1.0)
    
    #do it.  Want R*[n,e,d]
    n=localvec[0]
    e=localvec[1]
    d=localvec[2]
    
    n2=+r1*n+r2*e+r3*d
    e2=r4*n+r5*e+r6*d
    d2=r7*n+r8*e+r9*d
    
    #prep the answer for public consumption
    geovec=dim(n2,e2,d2)
    
    return geovec
           
#============================================
#geotolocal: Use this to take a vector from earth's center and recast it as
#a vector from earth's center
#
#Wants lon and lat of site on surface and lon and lat of the geocentric vector
#(where it pokes through surface)
#
#
def geotolocal(sitelon,sitelat,lon,lat):
    #set up the rotation matrix.
    #Cols are (local) NED axes in geocentric coords
    lam=radians(sitelon)   #spot on surface where we want local vector
    phi=radians(sitelat)
    r1=-sin(phi)*cos(lam)
    r2=-sin(lam)
    r3=-cos(phi)*cos(lam)
    r4=-sin(phi)*sin(lam)
    r5=cos(lam)
    r6=-cos(phi)*sin(lam)
    r7=cos(phi)
    r8=0.0
    r9=-sin(phi)
    
    #pole  to rec the site where we want local direction
    sitevec=ned(sitelon,sitelat,1.0)
    
    #pole to rec the geocentric direction
    geovec=ned(lon,lat,1.0)
    n=geovec[0]
    e=geovec[1]
    d=geovec[2]
    
    #Do it. Want R(transpose)*[n,e,d]
    n2=r1*n+r4*e+r7*d
    e2=r2*n+r5*e+r8*d
    d2=r3*n+r6*e+r9*d
    
    #rec to pole the answer
    localvec=dim(n2,e2,d2)
    
    return localvec
#=====================================================
#matmaker: Makes an array (list of row-lists) from a list
#Wants the list of elements and # of rows.
#Returns the array
def matmaker(elements,rows):
    #Makes an array (list of row-lists) from a list.
    #
    #Check that all elements are listed. Return empty list if not.
    if len(elements)%rows!=0:
        matout=[]
    else:
        matout=[]
        ncols=len(elements)//rows
        count=0
        for i in range(rows):
            thisrow=[]
            for j in range(ncols):
                thisrow.append(elements[count])
                count=count+1
            matout.append(thisrow)
    return matout                  

#======================
#takes and arrray (list of row-lists) and returns a list
def matunmaker(array):
    outlist=[]
    if len(array)>0:
        rows=len(array)
        cols=len(array[0])
        for i in range (rows):
           for j in range(cols):
               outlist.append(array[i][j])
    return(outlist)
#===================================================
#getrotmat: Get rotation matrix for rot about D,I
#Wants rot axis dec, rot axis inc, and rot
#Returns  list of 3 row lists
#
def getrotmat(axisdec,axisinc,rot):
    a=ned(axisdec,axisinc,1)
    n1=a[0]
    n2=a[1]
    n3=a[2]
    t=radians(rot)
    cost=cos(t)
    sint=sin(t)

    R11=cost+n1*n1*(1-cost)
    R12=n1*n2*(1-cost)-n3*sint
    R13=n1*n3*(1-cost)+n2*sint
    R21=n1*n2*(1-cost)+n3*sint
    R22=cost+n2*n2*(1-cost)
    R23=n2*n3*(1-cost)-n1*sint
    R31=n1*n3*(1-cost)-n2*sint
    R32=n2*n3*(1-cost)+n1*sint
    R33=cost+n3*n3*(1-cost)

    Rlist=[R11,R12,R13,R21,R22,R23,R31,R32,R33]
    R=matmaker(Rlist,3)
    return R
#===============================================
#dorotmat: Uses rot matrix to do an ARot.  But its CCW!
#Wants dec,inc,axisdec,axisinc, rot
#Returns [D,I,n,e,d]
def dorotmat(dec,inc,axisdec,axisinc,rot):
    b=ned(dec,inc,1)
    bcol=matmaker(b,3)
    R=getrotmat(axisdec,axisinc,rot)
    outvec=matmult(R,bcol)
    #get dec and inc, but slice out the M=1
    vec=dim(outvec[0][0],outvec[1][0],outvec[2][0])[:-1]
    vec.append(outvec[0][0])
    vec.append(outvec[1][0])
    vec.append(outvec[2][0])
    return vec
#==================================================
#The following is a set of routine that deal with
#matrices, which are lists of row-lists.
#Got these off the web Jan 2020
#
def transposeMatrix(m):
    return list(map(list,list(zip(*m))))

def getMatrixMinor(m,i,j):
    return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]

def getMatrixDeterminant(m):
    #base case for 2x2 matrix
    if len(m) == 2:
        return m[0][0]*m[1][1]-m[0][1]*m[1][0]

    determinant = 0
    for c in range(len(m)):
        determinant += ((-1)**c)*m[0][c]*getMatrixDeterminant(getMatrixMinor(m,0,c))
    return determinant

def getMatrixInverse(m):
    determinant = getMatrixDeterminant(m)
    #special case for 2x2 matrix:
    if len(m) == 2:
        return [[m[1][1]/determinant, -1*m[0][1]/determinant],
                [-1*m[1][0]/determinant, m[0][0]/determinant]]

    #find matrix of cofactors
    cofactors = []
    for r in range(len(m)):
        cofactorRow = []
        for c in range(len(m)):
            minor = getMatrixMinor(m,r,c)
            cofactorRow.append(((-1)**(r+c)) * getMatrixDeterminant(minor))
        cofactors.append(cofactorRow)
    cofactors = transposeMatrix(cofactors)
    for r in range(len(cofactors)):
        for c in range(len(cofactors)):
            cofactors[r][c] = cofactors[r][c]/determinant
    return cofactors

#elegant one-liner that does matrix multiplication
def matmult(m1, m2):
    return [
        [sum(x * y for x, y in zip(m1_r, m2_c)) for m2_c in zip(*m2)] for m1_r in m1
    ]
#===================================================
#eigen: calculate Eigenvectors and values for 3x3 matrix.
#
#13 Oct 2020   Added a test checking whether eigenvectors are
#orthogonal (within a degree).
#
#Will always return *something* now.
#If OK is false and there are numbers other than zero in the eigens,
#then probably you will have non-orthogonal eigenvectors.
#That will happen when 2 (or 3) of the eigenvalues are the same.
#
#Found a bug?  31 July 2020 Added "float" to calcs of a,b,c
#to avoid trouble with integers coming in.  Was not the problem!
#
#The input is a list of 6 numbers:
#[k11,k12,k13,k22,k23,k33]  upper half of symmetric 3x3 matrix
#Test data:
#           [2.3500  0.0239 -0.2400]
#           [        2.0800 -0.0014]
#           [                2.2200]
#
#           K1=2.5345  D=183.1 I=37.3
#           K2=2.0830  D= 82.1 I=14.2
#           K3=2.0325  D=335.1 I=49.2
#
#Output is a list of the 3 [eigen-dec,eigen-inc,eigen-val] lists
#Pythonized from the Pascal and lightly tested by Bogue 11-20-05
#8/7/2020  Checked against online eigenvector calulator - perfect match!
#
def eigen(inmat):
    from pmag import ang
    
    #--------------
    #this gives you +/- 1 millonth of av
    def tickle(av):
        from random import randint
        a=randint(0,1)
        if a>0:
            b=av*1E-6
        else:
            b=-av*1E-6
        return b    
    #-----------------
    #this is the old eigen, but now inside a try/except
    #eigout[3] is now a boolean to check for bombout
    def doit(m11,m12,m13,m22,m23,m33):
        ok=True
        eigout=[[0,0,0],[0,0,0],[0,0,0],ok]
        try:
            a=float(m11+m22+m33)
            b=float(m12*m12+m23*m23+m13*m13\
              -m11*m22-m22*m33-m11*m33)
            c=float( m11*m22*m33   \
              +2*m12*m23*m13 \
              -m12*m12*m33   \
              -m23*m23*m11   \
              -m13*m13*m22)

            pp=-b-(a*a/3)
            qq=-c-(a*b/3)-(2*a*a*a/27)
            rr=sqrt(-pp*pp*pp/27)
            theta=acos(-qq/(2*rr))
            CubeR=exp(log(rr)/3)   #log with 1 arg only is ln
            
            Eig1=2*CubeR*cos(theta/3)+a/3
            Eig2=2*CubeR*cos((theta+2*pi)/3)+a/3
            Eig3=2*CubeR*cos((theta+4*pi)/3)+a/3
            
            eiglist=[Eig1,Eig2,Eig3]
            eiglist.sort()     #sorts in place from low to hi
            eiglist.reverse()  #now hi to lo
                
            #eigout=list()  #empty list to append to in loop
            j=0  #0=first row, etc.
            while j<3:
                denom=m12*m23 + m13*(eiglist[j]-m22)
                N=1.0
                E=   (  (eiglist[j]-m11)*m23 + m13*m12   )/denom  
                D=   (  -m12*m12 + (eiglist[j]-m11) * (eiglist[j]-m22)  )/denom   
                eigdir=dim(N,E,D)  #so eigdir has dec,inc,mom
                if eigdir[1]<0:    #make eigenvectors point down
                    newI=-eigdir[1]
                    if eigdir[0]<180: newD=eigdir[0]+180
                    else:             newD=eigdir[0]-180
                    eigdir=[newD,newI]
                eigout[j]=[eigdir[0],eigdir[1],eiglist[j]]                
                j=j+1
        except:
            eigout[3]=False
                
        return(eigout)    
    #------------
    
    m11=inmat[0] #make it look like standard matrix 
    m12=inmat[1]
    m13=inmat[2]
    m22=inmat[3]
    m23=inmat[4]
    m33=inmat[5]

    
    eigout=doit(m11,m12,m13,m22,m23,m33)

    if not eigout[3]:
        print('It bombed!')

    itworked=eigout[3]

    if not itworked:
        #bombed!  Probably 2 eigenvalues equal (or 3)
        #Lets jiggle things a bit...
        av=(m11+m22+m33)/3.0
        m11=tickle(av)+m11
        m12=tickle(av)+m12
        m13=tickle(av)+m13
        m22=tickle(av)+m22
        m23=tickle(av)+m23
        m33=tickle(av)+m33
        #print (m11,m12,m13,m22,m23,m33)
        eigout=doit(m11,m12,m13,m22,m23,m33)

    itworked=eigout[3]

    if not itworked:
        print ('It bombed... even after tickling')    

    #No bombout, but did we get non-orthogonal eigen vector?

    a1=ang(eigout[0][0],eigout[0][1],eigout[1][0],eigout[1][1])
    a2=ang(eigout[0][0],eigout[0][1],eigout[2][0],eigout[2][1])
    a3=ang(eigout[1][0],eigout[1][1],eigout[2][0],eigout[2][1])

    if abs(a1-90)>1.0:
        eigout[3]=False
    if abs(a2-90)>1.0:
        eigout[3]=False    
    if abs(a3-90)>1.0:
        eigout[3]=False

    return eigout  #which is a list of 3 eigens [dec,inc,eigenval] + ok.

#====================================================================
#Function eigmat:  Calculate moment matrix from a list of
#directions in form [dec,inc,mom].
#Eigen dirs and vals of this matrix will be axes of cigar or pancake fit.
#Eigen vals give variance *along* selected axis.
#Input is the list of directions and the boolean free.
#if free=1=True, then cloud of vector endpoints will be shifted to be
#centered on its center of mass.
#If free=0=False, then fit to vector endpoints will be forced through origin.
#Pythonized from the Pascal and lightly tested by Bogue 11-20-05.
#Bombs if points in a perfect line or (?) plane.
def eigmat(dilist,free):
    n=e=d=nn=ee=dd=ne=nd=ed=0  #accumulators
    i=0
    num=len(dilist)
    #first pass thru list to get av values
    while i<num:
        thisdi=dilist[i]
        thisNED=ned(thisdi[0],thisdi[1],thisdi[2])
        n=n+thisNED[0]
        e=e+thisNED[1]
        d=d+thisNED[2]
        i=i+1
        
    #get the avs, pos of the center of mass    
    nav=n/num
    eav=e/num
    dav=d/num

    #second pass to fill up cos matrix, free or anchored
    i=0
    while i<num:
        thisdi=dilist[i]
        thisNED=ned(thisdi[0],thisdi[1],thisdi[2])
        thisN=thisNED[0]
        thisE=thisNED[1]
        thisD=thisNED[2]
        if free:   #shift cloud to center of mass
            thisN=thisN-nav
            thisE=thisE-eav
            thisD=thisD-dav
        nn=nn+thisN*thisN
        ee=ee+thisE*thisE
        dd=dd+thisD*thisD
        ne=ne+thisN*thisE
        nd=nd+thisN*thisD
        ed=ed+thisE*thisD
        i=i+1
    mommat=[nn,ne,nd,ee,ed,dd]  #ready to feed to eigen
    return mommat

#============================================
#eigen2x2
#Gives eigevector and eigenvalues of a 2x2 matrix.
#Input is a list of lists: [[11,12],[21,22]]
#That is:  [row1,row2]
#Returns [vec11,vec12,eigenvalue1],[vec21,vec22,eigenvalue2]
#The eigenvectors have unit length
#Shouts at you and returns 99999 if no real roots
def eigen2x2(A):
    a=1.0
    b=-(A[0][0]+A[1][1])
    c=A[0][0]*A[1][1]-A[1][0]*A[0][1]
    #print a,b,c
    if (b*b-4*a*c)>=0:
        root1=0.5*(-b+sqrt(b*b-4*c))
        root2=0.5*(-b-sqrt(b*b-4*c))
    else:
        print('NO REAL ROOTS!')
        root1=99999
        root2=99999
    #print root1, root2
    vec11=(A[1][1]-root1)/A[1][0]
    vec12=-((A[0][0]-root1)*vec11)/A[0][1]
    norm=(vec11*vec11+vec12*vec12)**0.5
    vec11=vec11/norm
    vec12=vec12/norm
    #print vec11,vec12
    #print sqrt(vec11*vec11+vec12*vec12)
    vec21=(A[1][1]-root2)/A[1][0]
    vec22=-((A[0][0]-root2)*vec21)/A[0][1]
    norm=(vec21*vec21+vec22*vec22)**0.5
    vec21=vec21/norm
    vec22=vec22/norm
    #print vec21,vec22
    #print sqrt(vec21*vec21+vec22*vec22)
    eigen1=[vec11,vec12,root1]
    eigen2=[vec21,vec22,root2]
    return [eigen1,eigen2]
    
#=============================================
#Use this calculate a Kirschvink line fit to paleomag data
#Input is a list if [D,I,M,A] elements and boolean "free"
#Note: A ignored   
#If free is False, anchored line fit
#Returns the fit in list [fitD,fiti,MAD,pathindex]
#New 8/2020 "pathindex"
#Pathindex one when points form a line in sequence.
#Pathindex<1 if path wanders or backtracks.
#Lightly tested March 2011 by swb

def linefit(inlist,free):
    #from math import sqrt,degrees,atan2
    thismat=eigmat(inlist,free)
    thiseiglist=eigen(thismat)
    '''
    i=0
    while i<3:
        print(thiseiglist[i])
        i=i+1
    '''
    fitD=thiseiglist[0][0]
    fitI=thiseiglist[0][1]

    #Check if we have first ev pointing the right way.
    #Take first vector minus last and see if it is
    #within 90 deg of eig vector. If not, invert.

    #slice out first 3 - D,I,M
    vec1=inlist[0][0:3]
    nedlist1=ned(vec1[0],vec1[1],vec1[2])
    #print vec1
    vec2=inlist[len(inlist)-1][0:3]
    nedlist2=ned(vec2[0],vec2[1],vec2[2])
    #print vec2
    #print nedlist2
    diffned=[nedlist1[0]-nedlist2[0],nedlist1[1]-nedlist2[1],nedlist1[2]-nedlist2[2]]
    diffdir=dim(diffned[0],diffned[1],diffned[2])
    #print diffdir
    #now get nedlist for eigvector
    nedlistev=ned(fitD,fitI,1)
    #Get angle between em
    dotlist=dot(diffned,nedlistev)
    ang=dotlist[1]
    #Do inversion check
    if ang>90.0:
        fitD=fitD-180.0
        if fitD<0:
            fitD=fitD+360.0
        fitI=-fitI
    #print fitD,fitI
    #get max angle of deviation
        
    #min eig values could be negative if really 0.
    thiseiglist[1][2]=abs(thiseiglist[1][2])   
    thiseiglist[2][2]=abs(thiseiglist[2][2])   

    #this gives same answer as Craig Jones program
    mad=degrees(atan2(sqrt(thiseiglist[1][2] + thiseiglist[2][2]),sqrt(thiseiglist[0][2])))

    #Get the ratio of vecdiff(last-first)/arithmetic sum of mom diffs.   Will be close to 1 if
    #points are collinear and path is seqential; will be <1 if the points backtrack a lot - path wanders.
    #Possible to get low MAD but have this ratio<<1.
    #Already got first-last above (diffdir)
    arithsum=0
    for i in range(1,len(inlist)):
        thisone=inlist[i][0:3]
        thatone=inlist[i-1][0:3]
        thisned=nedlist(thisone)
        thatned=nedlist(thatone)
        diff=dimlist([thisned[0]-thatned[0],thisned[1]-thatned[1],thisned[2]-thatned[2]])
        arithsum=arithsum+diff[2]
    pathindex= diffdir[2]/arithsum          
    
    return [fitD,fitI,mad,pathindex]

#==============================================
#Find plane or great circle that best fits a set of points
#
#Takes a list of [D,I] items
#
#Returns a list of [pole dec, pole inc, eigenvalue,n]
def planefit(inlist):
    n=len(inlist)
    #turn each [D,I] into [D,I,M]
    veclist=list()
    for i in range(n):
        thisitem=inlist[i]
        thisitem.append(1.0)
        veclist.append(thisitem)
    mommat=eigmat(veclist,False)
    eigstuff=eigen(mommat)
    #get the one with lowest eigenvalue
    mineig=min(eigstuff[2])
    for i in range (3):
        if eigstuff[i][2]==mineig:
            pdec=eigstuff[i][0]
            pinc=eigstuff[i][1]
    return [pdec,pinc,mineig,n]

#==============================================
#Find point on great circle closest to a direction
#
#Takes(dec,inc,poledec,poleinc) and returns [dec2,inc2,ang]
#Note: ang is angle from dir to gc
#

def closest(dec,inc,poledec,poleinc):
    #get the direction into rect
    ned1=ned(dec,inc,1.0)
    #get the pole to gc into recgt
    ned2=ned(poledec,poleinc,1.0)
    #cross the two
    ned3=crossunit(ned1,ned2)
    #now get the point on the gc closest to ned1
    ned4=crossunit(ned2,ned3)
    dim4=dimlist(ned4)
    #print dim4
    xpdec=dim4[0]  #xp = crossing point
    xpinc=dim4[1]
    dot1=dot(ned1,ned4)
    ang=dot1[1]
    return [xpdec,xpinc,ang]

#====================================================
#Orthogonal regression  (aka Deming regression)
#
#Assumes equal errors on both x and y
#Find line that minimizes sum of orthgonal distances to data
#
#Wants two lists of same length: xvals and yvals
#
#Returns slope and y-intercept
#
#added 6 Jan 2017
#
def orthoreg(xvals,yvals):
    #find means
    i=0
    xsum=0
    ysum=0
    n=len(xvals)
    #print n
    for i in range(n):
        xsum=xsum+xvals[i]
        ysum=ysum+yvals[i]
    N=n*1.0  #make it real before dividing     
    xm=xsum/N
    ym=ysum/N
    #print xm, ym
    sxx=0
    sxy=0
    syy=0
    for i in range (n):
        sxx=sxx+(xvals[i]-xm)*(xvals[i]-xm)
        sxy=sxy+(xvals[i]-xm)*(yvals[i]-ym)
        syy=syy+(yvals[i]-ym)*(yvals[i]-ym)
    sxx=sxx/(N-1)
    sxy=sxy/(N-1)
    syy=syy/(N-1)
    #print sxx,sxy,syy
    m=syy-sxx+((syy-sxx)*(syy-sxx)+4*sxy*sxy)**(0.5)
    m=m/(2*sxy)
    b=ym-m*xm
    outlist=(m,b)
    return outlist

#====================================================
#All regressions(least squares on y, x, and orthogonal) 
#
#Wants two lists of same length: xvals and yvals
#
#Returns 3 slopes and, someday, 3 y-intercepts
#And, someday, uncertainties on all 3 slopes.
#
#added 6 Jan 2017
#
def allreg(xvals,yvals):
    #find means
    i=0
    xsum=0
    ysum=0
    xxsum=0
    yysum=0
    xysum=0
    n=len(xvals)
    #print n
    for i in range(n):
        xsum=xsum+xvals[i]
        ysum=ysum+yvals[i]
        xxsum=xxsum+xvals[i]**2
        yysum=yysum+yvals[i]**2
        xysum=xysum+xvals[i]*yvals[i]
    N=n*1.0  #make it real before dividing     
    xm=xsum/N
    ym=ysum/N
    #print xm, ym
    sxx=0
    sxy=0
    syy=0
    for i in range (n):
        sxx=sxx+(xvals[i]-xm)*(xvals[i]-xm)
        sxy=sxy+(xvals[i]-xm)*(yvals[i]-ym)
        syy=syy+(yvals[i]-ym)*(yvals[i]-ym)
    sxx=sxx/(N-1)
    sxy=sxy/(N-1)
    syy=syy/(N-1)
    #Get params for orthogonal
    #print sxx,sxy,syy
    m=syy-sxx+((syy-sxx)*(syy-sxx)+4*sxy*sxy)**(0.5)
    m=m/(2*sxy)
    b=ym-m*xm
    #Get params for least squares on y
    my=(N*xysum-(xsum*ysum))/(N*xxsum-xsum**2)
    #Get params for least squares on x
    mx=(N*xysum-(xsum*ysum))/(N*yysum-ysum**2)
    mx=1/mx
    outlist=(m,my,mx)
    return outlist

#==================================================================
#ListStats(inlist)

#Returns mean, sd, and sd/mean and n of list
#
#
def ListStats(inlist):
    n=len(inlist)
    sum=0
    for i in range(n):
        sum=sum+inlist[i]
    mean=sum/(1.0*n)    
    sumsq=0
    for i in range(n):
        sumsq=sumsq+(mean-inlist[i])*(mean-inlist[i])
    sd=sqrt(sumsq/(n-1.0))
    sdm=sd/mean
    outlist=[mean,sd,sdm,n]
    return outlist

#==========================================================
#thellfit(xvals,yvals)
#
#Does the Coe et al. 1983 line fit that weights x and y residuals by
#their respective variances. Gives a slope and standard error of the slope.
#
#Note well: can't tell a positive slope from a negative slope!
#So check that with pmag.orthoreg
#
#So this gives a negative sb/b if b is neg.  Does any other program care?
#If not, would be good to change sign if b is neg.
#
#Wants two lists: xvals and yvals.
#Returns a list with b,sbsb/b,N
#
def thellfit(xvals,yvals):
    
    #check if slope is + or -
    test=orthoreg(xvals,yvals)
    if test[0]<0.0:
        negslope=True
    else:
        negslope=False
        
    n=len(xvals)
    sumx=0
    sumy=0
    for i in range(n):
        sumx=sumx+xvals[i]
        sumy=sumy+yvals[i]
    meanx=sumx/(1.0*n)
    meany=sumy/(1.0*n)
    sumxx=0
    sumyy=0
    sumxy=0
    for i in range(n):
        sumxx=sumxx+(xvals[i]-meanx)*(xvals[i]-meanx)
        sumyy=sumyy+(yvals[i]-meany)*(yvals[i]-meany)
        sumxy=sumxy+(xvals[i]-meanx)*(yvals[i]-meany)

    if negslope:    
        b=-sqrt(sumyy/sumxx)
    else:
        b=sqrt(sumyy/sumxx)
        
    sb2=(2*sumyy-2*b*sumxy)/((n-2.0)*sumxx)
    sb=sqrt(sb2)

    #might as well get the y-intercept too
    yint=meany-(b*meanx)
    outlist=[b,sb,sb/b,n,yint]
    return outlist


#==================================================================
#GetLev(min,max,n)  Finds n-1 contour intervals centered on zero.
#                   Gives n contour lines, half above zero.
#                   Make sure n is odd!
#
#Returns a list ready for dislin.contur
#
#n=number of contour lines - make it odd to get zero contour
#assuming we want contours centered on zero, like for radial field
#returns a list with the contour values

def GetLev(min,max,n):
    levlist=list()  #put answers here
    
    if max>=-min:
        absmax=max
    else:
        absmax=-min
        
    #get range up where we can see it
        
    up=True

    countdn=0
    levrange=2.0*absmax
    while levrange<10:
        up=False
        levrange=10.0*levrange
        countdn=countdn+1
        
    #so now between 10 and 99.9999 or >100

    countup=0
    while levrange>100:
        levrange=levrange/10.0
        countup=countup+1

    #so now between  10 and 99.999

    #get the interval   
    step=0.1
    while step*(n-1)<levrange:
        step=step+0.1
        
    #rescale them to levrange
    if up:
        step=step*10.0**countup
    else:
        step=step/10.0**countdn
        
    #get the levs.  Int division n\2 truncates
    startlev=0.0-((n/2)*step)
    for i in range (n):
        item=startlev+i*step
        levlist.append(item)
        #print item        

    return levlist
    #print 'countdn,countup,up',countdn,countup, up     
       
   
