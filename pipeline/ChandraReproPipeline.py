#!/usr/bin/python

# ==========================================
# ==========================================
# NAME:
#   ChandraReproPipeline
# 
# PURPOSE: 
#   To download, reprocess, and perform basic data extraction
#   for Chandra SgrA* / Magnetar observations
#
# INPUTS:
#   Pipeline Parameter file string (including path)
#   descriptions of the variables are found within the parameter file
# 
# OUTPUTS:
#   Creates xxx
#
# EXAMPLES:
#   >> import ChandraReproPipeline as cpipe
#   >> cpipe.Run("/Full/Path/To/File/cpipe_parfile_15043.py")
#        or
#   $ python ~/scripts/ChandraReproPipeline.py par_15043.py 
#
# NOTES:
#   CIAO needs to be initialized PRIOR to starting python
# ==========================================
# ==========================================

import numpy as np
from astropy.io import fits
from astropy import wcs
import sys
import csv
from ChPipe_parfile import *

import ciao_contrib.runtool as rt
import glob
import sys
import imp
from subprocess import call
import os
import datetime
import ciao_contrib.logger_wrapper as lw
from ciao_contrib.cda.data import download_chandra_obsids
import shutil
import pychips as pc
import pychips.hlui as hlui
#from astropy.io import fits
#import matplotlib.pyplot as plt

wav_regions=np.zeros((1,1))
filtered_rad_reg=np.zeros((1,2))



# ==========================================
def Run(argv):

    print( "================================================ ")
    print( "    Beginning Chandra reprocessing pipeline")
    print( "================================================ ")
    
    # --------- import parameters from parfile
    if type(argv) is str:
        parfile = argv
    if type(argv) is list:
        parfile = argv[0]
    
    parfilename = str(parfile.split("/")[-1:][0])
    par = imp.load_source(parfilename, parfile)
        
    print( "The following parameter file was successfully uploaded:")
    print( "     "+parfile)
    parholder = {key: value for key, value in par.__dict__.items() if not (key.startswith('__') or key.startswith('_'))}
    for key, value in parholder.items():
        print( str(key)+" = "+str(value))
    
    # -------  these are some parameters that do not go into the parFile
    halfwidth = 90    # half width of box used to create cut evt2 file, in pix
    boxwidth = 88     # in arcsec
    
    # --------- beginning parameter processing
    # create erange
    eRange = par.eRange
    emin = int(eRange[0] * 1000)
    emax = int(eRange[1] * 1000)
    if emin % 10 == 0:
        etag = str(int(eRange[0])) + '-' + str(int(eRange[1])) + 'keV'
    else:
        etag = str(eRange[0]) + '-' + str(eRange[1]) + 'keV'
    # create plistFile
    psListFile = par.plistFile
    # create tags for the extracted lightcurves
    sradtag = str(par.sgrArad)+"asec"
    mradtag = str(par.magrad)+"asec"

    # start download process if required
    if par.downloadBool=='yes' or par.downloadBool=='YES' or par.downloadBool=="Yes":
        print( "-----------------")
        print( "beginning download process")
        print( "-----------------")
        downloadData(par.obsID_list, par.downloadPath)
                             
    # creates the repro directory as repro_Nov17
    tstamp= datetime.datetime.now()
    reproDir = "repro_"+tstamp.strftime("%b")+tstamp.strftime("%y")
    
    for obsID in par.obsID_list:
        print( "-----------------")
        print( "beginning analysis for "+str(obsID))
        print( "-----------------")
        
        # first process the requred full paths
        reproPath = par.reproPathBase + "/" + reproDir + "/" + str(obsID) + "/"
        dlPathFull = par.downloadPath+"/"+str(obsID)
        dlPrimPath = par.downloadPath+"/"+str(obsID)+"/primary/"

        print( "** creating directory for reprocessing data")
        if not os.path.exists(reproPath):
                os.makedirs(reproPath)
        print( "   "+reproPath)

        if par.reproBool=='yes' or par.reproBool=='YES' or par.reproBool=="Yes":
            print( "** reprocessing data")
            rt.chandra_repro(indir=dlPathFull, outdir=reproPath, root='acisf'+str(obsID), badpixel='yes',
                             process_events='yes', set_ardlib='yes', pix_adj='edser', check_vf_pha='no',
                             cleanup='yes', destreak='yes', clobbe='yes')
            print( "   finished reprocessing")

        print( "** finding the event, asol, mask, and badpixel files in reprecessed dataset")
        evt2File, asol1,  asol1lis, msk1, bpix  = findFiles(reproPath)   #contain the full paths and file names
        hkfiles = [asol1,  asol1lis, msk1, bpix]

        print( "** setting the badpix file ")
        setbpixFile(reproPath, bpix)

        print( "** converting SrgrA, magnetar, and background  ra/dec into x/y coordinates")
        [sgrApix, magpix, backpix] = getPixelCoords(reproPath, evt2File, [par.sgrApos, par.magpos, par.bkgpos], asol1)
        
        print( "** creating the energy- and spatial-cut event file for wavdetect")
        squareRegFile = makeRegFiles(reproPath, str(obsID)+"_wvdet_square.reg", sgrApix, [boxwidth, boxwidth], "box")
        squareEvtFile = dmcopyEvtFile(evt2File, reproPath, squareRegFile, obsID, "_square_2-8keV", erange=[emin, emax], ftype='evt')
        
        print( '** starting getPSFmap')
        psfmap = getPSFmap(squareEvtFile)

        print( "** starting wavdetect call")
        wavdetsrcFile = wvdetectCall(squareEvtFile, psfmap)
        
        print( "** creating regions and ds9 image before wcs corrections")
        sgrAregFile = makeRegFiles(reproPath, str(obsID)+"_sgrA_orig.reg", sgrApix, par.sgrArad, "circ")
        magregFile = makeRegFiles(reproPath, str(obsID)+"_magnetar_orig.reg", magpix, par.magrad, "circ")
        getds9image([sgrAregFile, magregFile], reproPath+ str(obsID)+"_repro_2-8keV_prewcs.jpg", squareEvtFile)

        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #adding alex's code 

        #extrcting the CRVAL values from the image files


        #getting the source list from wavdetect
        source_list=wavdetect_source(wavdetsrcFile)[0]
        print(source_list)
        all_sources=wavdetect_source(wavdetsrcFile)[1]

        #executing translation code 
        trans=translation(source_list,all_sources)

        #loading the regions from the ascii radio file and defining variables that are used in the function

        reg=load_regions()
        
        #filtering the radio region data file and assigning them to variables, filtering the wavdetect regions as well
        #declaring variables that will be used in the functions


        global filtered_rad_reg
        global wav_regions
        
        filtered_radio_data=filter_regions(all_sources,reg)    
        wav_regions=np.ndarray.tolist(np.delete(wav_regions,0,0))

        filtered_wav_data=np.column_stack([np.take(all_sources[...,0],wav_regions),np.take(all_sources[...,1],wav_regions)])
        print("On the left is the radio coordinates, and on the right are the wavedetect x-ray coordinates.")
        print(np.column_stack([filtered_radio_data,filtered_wav_data]))

        #Taking the offset of the distance of all sources

        mean_offset=all_translation(filtered_radio_data,filtered_wav_data)

        #saving the region files into the current directory and in the specified directory
        
        save_reg(filtered_radio_data,obsID,reproPath)

        #copying the files in order to manipulate the coordinate system 
        evt2wcsFile,asolwcs=copy_files(evt2File,asol1,reproPath,obsID)
        
        #modifying the acis event file 
        
        full_evt_translation=TCRVL(evt2wcsFile,mean_offset,reproPath)
        
        
        #modifying the asol file

        asol1_translation=asol1_update(asolwcs,mean_offset,reproPath)

        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        """
        commenting out wcs_routines function
        print( "** performing wcs corrections")
        evt2wcsFile, asolwcs = wcsRoutines(obsID, squareEvtFile, evt2File, psListFile, wavdetsrcFile, reproPath, asol1)
        print( "   "+evt2wcsFile)
        """
        hkfiles[0] = asolwcs
        
        print( "** creating regions and ds9 image after wcs corrections")
        [sgrApixwcs, magpixwcs, backpixwcs] = getPixelCoords(reproPath, evt2wcsFile, [par.sgrApos, par.magpos, par.bkgpos], asol1)
        sqReg_wcsFile = makeRegFiles(reproPath, str(obsID)+"_wvdet_square_wcs.reg", sgrApixwcs, [boxwidth, boxwidth], "box")
        sgrAreg_wcsFile = makeRegFiles(reproPath, str(obsID)+"_sgrA_wcs.reg", sgrApixwcs, par.sgrArad, "circ")
        magreg_wcsFile = makeRegFiles(reproPath, str(obsID)+"_magnetar_wcs.reg", magpixwcs, par.magrad, "circ")
        backreg_wcsFile = makeRegFiles(reproPath, str(obsID)+"_bkg_wcs.reg", magpixwcs, par.bkgrad, "ann")
        cutEvtwcsFile = dmcopyEvtFile(evt2wcsFile, reproPath, sqReg_wcsFile , obsID, "_square_2-8keV_wcs", erange=[emin, emax], ftype='evt')
        getds9image([sgrAreg_wcsFile, magreg_wcsFile,backreg_wcsFile], reproPath+str(obsID)+"_repro_2-8keV_postwcs.jpg", cutEvtwcsFile)


        # -----------
        print( "** beginning event file region extractions")
        sgrEcutwcsEvt = dmcopyEvtFile(evt2wcsFile, reproPath, sgrAreg_wcsFile , obsID, "_SgrA_2-8keV_wcs_evt", erange=[emin,emax], ftype='evt')
        sgrwcsEvt     = dmcopyEvtFile(evt2wcsFile, reproPath, sgrAreg_wcsFile , obsID, "_SgrA_wcs_evt", ftype='evt')
        magEcutwcsEvt = dmcopyEvtFile(evt2wcsFile, reproPath, magreg_wcsFile , obsID, "_magnetar_2-8keV_wcs_evt", erange=[emin,emax], ftype='evt')
        magwcsEvt     = dmcopyEvtFile(evt2wcsFile, reproPath, magreg_wcsFile , obsID, "_magnetar_wcs_evt", ftype='evt')
        
        print( "** beginning the light curve extractions for 2-8 keV")
        # the outputs here are lc file names"
        sgr_LCnobkg, sgr_LCbkg = makeLightCurves(sgrAreg_wcsFile, backreg_wcsFile, sradtag, "sgra", sgrEcutwcsEvt, obsID, reproPath,
                                                 emin, emax, etag, par.lcbin, par.ccd_id)
        mag_LCnobkg, mag_LCbkg = makeLightCurves(magreg_wcsFile, backreg_wcsFile,  mradtag, "magnetar", magEcutwcsEvt, obsID, reproPath,
                                                 emin, emax, etag, par.lcbin, par.ccd_id)
        print( "** plotting lightcurves")
        plotCiaoLcurve(sgr_LCnobkg, sgr_LCbkg, "sgrA", par.lcbin, obsID, reproPath)
        plotCiaoLcurve(mag_LCnobkg, mag_LCbkg, "magnetar", par.lcbin, obsID, reproPath)
        
        
        print( "** creating spectral file"        )
        sgrSpecFile = getSpectrum(evt2wcsFile, sgrAreg_wcsFile, "sgra", sradtag, hkfiles, obsID, reproPath, bkgReg=backreg_wcsFile)
        magSpecFile = getSpectrum(evt2wcsFile, magreg_wcsFile, "magnetar", mradtag, hkfiles, obsID, reproPath, bkgReg=backreg_wcsFile)
        
        
# ======================================
# ======================================
def downloadData(idList, downloadPath):
        
    currentDir = os.getcwd()
    os.chdir(downloadPath)
    print( "downloading data to the download dir")
    print( "    --> "+downloadPath)

    # this controls verbosity of download
    lw.initialize_logger('download', verbose=1)
    # Downloads this to the download (local) directory
    for idNumb in idList:
        download_chandra_obsids([idNumb])
        # if needed, unzip the files
        print( "unzipping")
        os.system("gunzip -r "+downloadPath+"/"+str(idNumb) )
        
    # change back to the original directory    
    os.chdir(currentDir)

    
# ======================================
def findFiles(dataPath):

    foundfiles = glob.glob(dataPath + '*repro_evt2.fits')
    if len(foundfiles) == 1:
        evt2file = foundfiles[0]
    else:
        print( 'error: either none or many evt2 files have been found')
        quit()
    foundfiles = glob.glob(dataPath + '*asol1.fits')
    if len(foundfiles) == 1:
        asol1file = foundfiles[0]
    else:
        print( 'error: either none or many asol1.fits files have been found')
        quit()
    foundfiles = glob.glob(dataPath + '*asol1.lis')
    if len(foundfiles) == 1:
        asol1lis = foundfiles[0]
    else:
        print( 'error: either none or many asol1.lis files have been found')
        quit()
    foundfiles = glob.glob(dataPath + '*msk1.fits')
    if len(foundfiles) == 1:
        msk1file = foundfiles[0]
    else:
        print( 'error: either none or many msk files have been found')
        quit()
    foundfiles = glob.glob(dataPath + '*repro_bpix1.fits')
    if len(foundfiles) == 1:
        bpixfile = foundfiles[0]
    else:
        print( 'error: either none or many repro_bpix1 files have been found')
        quit()

    return evt2file, asol1file, asol1lis, msk1file, bpixfile

# ======================================
def setbpixFile(dataPath, bpixFile):
    rt.acis_set_ardlib.punlearn()
    rt.acis_set_ardlib(bpixFile)
    return

# ======================================        
def getPixelCoords(dataPath, evt2, posArray, asol1): #sPos, mPos, bPos):            
    
    outArray = []
    for position in posArray:
        rt.dmcoords(evt2, asolfile=asol1, ra=position[0], dec=position[1], option='cel', verbose=1)
        pix = [rt.dmcoords.x, rt.dmcoords.y]
        print( "   pix     "+str(pix[0])+", "+str(pix[1]))
        outArray.append(pix)

    return outArray
    '''
    rt.dmcoords(evt2, asolfile=asol1, ra=sPos[0], dec=sPos[1], option='cel', verbose=1)
    sPix = [rt.dmcoords.x, rt.dmcoords.y]
    print( "   sgrA pix     "+str(sPix[0])+", "+str(sPix[1]))
    rt.dmcoords(evt2, asolfile=asol1, ra=mPos[0], dec=mPos[1], option='cel', verbose=1)
    mPix = [rt.dmcoords.x, rt.dmcoords.y]
    print( "   magnetar pix "+str(mPix[0])+", "+str(mPix[1]))
    rt.dmcoords(evt2, asolfile=asol1, ra=bPos[0], dec=bPos[1], option='cel', verbose=1)
    bPix = [rt.dmcoords.x, rt.dmcoords.y]
    print( "   bkg pix      "+str(bPix[0])+", "+str(bPix[1]))
    '''
    return sPix, mPix, bPix            

# ======================================
def makeRegFiles (dataPath, regFile, coords, radius, shape ):

    fullregFile = dataPath+regFile
    f = open(fullregFile, 'w')
    f.write("# Region file format: CIAO version 1.0 \n")
    if shape=="circ":
        rad = str(radius)+"\""
        f.write("circle(%s,%s,%s)" %(coords[0],coords[1],rad) )
    if shape=="ann":
        rad1 = str( radius[0])+"\""
        rad2 = str( radius[1])+"\""
        f.write("annulus(%s,%s,%s,%s)" %(coords[0],coords[1],rad1,rad2))
    if shape=="box":
        width = str( radius[0])+"\""
        height= str( radius[1])+"\""
        f.write("box(%s,%s,%s,%s,0)"  %(coords[0],coords[1],width,height))
    f.close()
    print( "      created "+regFile)
    return fullregFile


# ======================================
def getds9image( regArray, imgName, evtFile):

    # this creates a command to pass to the os.
    # creates a jpg image that overlays the regions on the evt File image
    ds9Routine = "ds9 "+evtFile+" "
    # append the regions if they exist
    if regArray!="none" or regArray!="NONE":
        for reg in regArray:
            ds9Routine = ds9Routine+" -regions load "+reg
    # append the image formatting parameters 
    ds9Routine = ds9Routine + " -zoom 8 -scale log -cmap b"
    # append the image name "
    ds9Routine = ds9Routine + " -saveimage jpeg " + imgName
    ds9Routine = ds9Routine + " -exit"
    # call ds9 to the command line
    os.system(ds9Routine)
    return

# ======================================
def dmcopyEvtFile(evt2, dataPath, region, obsID, tag, erange=[], ftype='evt'):

    inEvtFile = evt2+"[EVENTS]"+"[sky=region("+region+")]"
    if len(erange)>0:
        inEvtFile = inEvtFile+"[energy="+str(erange[0])+":"+str(erange[1])+"]"
    opt = 'all'
    if ftype=='img':
        opt='image'
    outEvtName = str(obsID)+tag+'.fits'
    outEvtFile = dataPath+outEvtName
    rt.dmcopy(infile=inEvtFile, outfile=outEvtFile, option=opt, clobber='yes') 
        
    print( "      created "+outEvtName)
    return outEvtFile
    
# ======================================
def getPSFmap(cutevtFile):
    
    psfmap = cutevtFile[:-5] + '_psfmap.fits'
    rt.mkpsfmap(cutevtFile, outfile=psfmap, energy=3.8, ecf=0.393, clobber='yes')
    return psfmap

# ====================================== 
def wvdetectCall(cutevtFile, psfmap):

    outFile = cutevtFile[:-5] + '_src.fits'
    scellFile = cutevtFile[:-5] + '_scell.fits'
    imageFile = cutevtFile[:-5] + '_img.fits'
    nbkgFile = cutevtFile[:-5] + '_nbkg.fits'
    regFile = cutevtFile[:-5] + '_src.reg'
    scales = '1.0 2.0 4.0 8.0 16.0'
    sigThresh = 1e-06
    rt.wavdetect(cutevtFile, psffile=psfmap, outfile=outFile,
                 scellfile=scellFile, imagefile=imageFile,
                 defnbkgfile=nbkgFile, regfile=regFile,
                 scales=scales, sigthresh=sigThresh, clobber='yes')
    return outFile

# ++++++++++++++++++++++++++++++++++++++

#change the file names, convert sexa coords to deg 
#inserting Alex's code for manual change in coordinate system


#=================================================
def deg_to_arc(x):
    y=x*3600
    return y 

#=================================================
#locates the position of the image tangent point
def CRVAL(filename):

    #loads the fits hdulist 
    hdulist=fits.open(filename, mode='update')
    
     # Parse the WCS keywords in the primary HDU
    w = wcs.WCS(hdulist[0].header)

    # Print out the RA and DEC position of the tangent point
    
    CRVAL1= w.wcs.crval[0]
    CRVAL2 =w.wcs.crval[1]

    return CRVAL1,CRVAL2

#================================================
#this function anaylyses the source file and takes out the three brightest sources detected in wav detect and places them in arrays.
def wavdetect_source(src_file):

        #loads the fits hdulist 
    hdulist=fits.open(src_file)

    wav=np.array((hdulist[1].data['RA'],hdulist[1].data['DEC']))

    net_counts=(hdulist[1].data['net_counts'])

    #to return where in the array the maximum value is, use argmax, to return more than one max value, use argsort.

    max= np.array(net_counts.argsort()[-3:][::-1])

    #this array contains the information of the coordinates of the 3 brightest sources in the given file
    
    source_list=np.column_stack([wav[0][max],wav[1][max],net_counts[max]])
    all_sources=np.column_stack([wav[0][:],wav[1][:],net_counts[:]])
    
    return source_list,all_sources
#In source_list, the RA coordinates are in the [:,0]
#In source_list, the DEC coordinates are in the [:,1]
#In source_list, the net_counts are in the [:,2]

#===============================================
#the translation of the wcs_coordinate system
def translation(src,allsources):
    
    while_var=False
    counts=0
    def in_circle(center_x, center_y, x, y,radius):

        
        square_dist = (center_x - x) ** 2 + (center_y - y) ** 2
        if square_dist <= radius ** 2:

            sag_x_ra=x
            sag_x_dec=y
            while_var=True
            
        else:
            sag_x_ra="Error"
            sag_x_dec="Did not detect source close to radio position of Sag A*"
            while_var=False
        return while_var,sag_x_ra,sag_x_dec

#    from ChPipe_parfile.py import sgrApos,magpos

    while not while_var and counts+1<len(allsources):
        counts=counts+1
        while_var=in_circle(sexa_to_deg(sgrApos[0]),sexa_to_deg(sgrApos[1]),allsources[counts][0],allsources[counts][1],2/3600)[0]
        

    sag_x_ra=in_circle(sexa_to_deg(sgrApos[0]),sexa_to_deg(sgrApos[1]),allsources[counts][0],allsources[counts][1],2/3600)[1]
    sag_x_dec=in_circle(sexa_to_deg(sgrApos[0]),sexa_to_deg(sgrApos[1]),allsources[counts][0],allsources[counts][1],2/3600)[2]
    print("Sag A* wavdetect output")
    print(sag_x_ra,sag_x_dec)


    ra_mag_offset=sexa_to_deg(magpos[0])-src[0][0]
    dec_mag_offset=sexa_to_deg(magpos[1])-src[0][1]

    ra_sag_offset=sexa_to_deg(sgrApos[0])-sag_x_ra
    dec_sag_offset=sexa_to_deg(sgrApos[1])-sag_x_dec

    average=np.array([(ra_mag_offset+ra_sag_offset)/2,(dec_mag_offset+dec_sag_offset)/2])
    print("Mag RA,DEC offset")
    print(deg_to_arc(ra_mag_offset))
    print(deg_to_arc(dec_mag_offset))
    print("SagA* RA,DEC offset")
    print(deg_to_arc(ra_sag_offset))
    print(deg_to_arc(dec_sag_offset))
    print("This is average offsest for both sources.")
    print(deg_to_arc(average))

    return average

#======================================
#changing the image tangent point
def getDmhedit(filename,CRVAL_num,value):

    infile=filename
    filelist='none'
    operation='add'
    key=CRVAL_num
    value=value
    datatype= 'indef'
    clobber='yes'
    
   
    
    rt.dmhedit.punlearn()
    rt.dmhedit(infile=infile, filelist=filelist, operation=operation, key=key, value=value, datatype=datatype)

#========================================================
#loading the ascii file into a region file
def load_regions():
#    from ChPipe_parfile.py import supplementalPath
    os.chdir(supplementalPath)


    with open(plistFile, 'rt') as f:
        reader = csv.reader(f, delimiter=' ', skipinitialspace=True)

        lineData = list()
        cols = next(reader)
        print(cols)

        for col in cols:
            # Create a list in lineData for each column of data.
            lineData.append(list())


        for line in reader:
            for i in range(0, len(lineData)):
                # Copy the data from the line into the correct columns.
                lineData[i].append(line[i])
            
        data = dict()

        for i in range(0, len(cols)):
            # Create each key in the dict with the data in its column.
            data[cols[i]] = lineData[i]

    # adding words to a list 
    reg=np.array(lineData)
    reg=reg.astype(np.float)

    reg_ra=reg[:][0]
    reg_dec= reg[:][1]
    return reg

"""
===================================================================
This function takes the outputs from the wavdetect source file and
filters the radio region catalogue file in the input. 

===================================================================
"""

#=================================================
#filtering the radio data with the x-ray positions
def filter_regions(source_list,regions):

    #============================
    #detecting points in circle around radio regions
    def in_circle(center_x, center_y, x, y,radius):


        global filtered_rad_reg
        global wav_regions
        square_dist = (center_x - x) ** 2 + (center_y - y) ** 2
        if square_dist <= radius ** 2:

            filtered_rad_reg=np.concatenate((filtered_rad_reg,np.array([[x,y]])))

            wav_regions=np.array(np.append(wav_regions,i))
          
        return filtered_rad_reg

    for i in range (0,len(source_list[...,0])):
            
       for j in range (0,len(regions[:][0])):
            coords=[source_list[i][0],source_list[i][1],regions[0][j],regions[1][j],1/3600]
            filtered_rad_reg=in_circle(*coords)
            
    filtered_rad_reg=np.delete(filtered_rad_reg,0,0)        
    return filtered_rad_reg

#==============================================
#taking the avergage translation of all sources
def all_translation(radio_data,xray_data):
    mean_offset=np.array([np.mean(radio_data[...,0]-xray_data[...,0]),np.mean(radio_data[...,1]-xray_data[...,1])])

    #if i install pandas, then the following command will be executable to make nice tables
    #pd.DataFrame([(deg_to_arc(filtered_radio_data[...,0]-filtered_wav_data[...,0])),(deg_to_arc(filtered_radio_data[...,1]-filtered_wav_data[...,1]))],["RA offset","DEC offset"])

    print("RA offset for each source")
    print(deg_to_arc(radio_data[...,0]-xray_data[...,0]))
    print("DEC offset for each source")
    print(deg_to_arc(radio_data[...,1]-xray_data[...,1]))
    print("The average offset is: ")
    print(deg_to_arc(mean_offset))
    return mean_offset

#===========================================
#saving the regions in the current PATH
def save_reg(radio_data,obsID,Path):
    global filtered_radio_data
    def converter(list1,list2):
        return "circle("+str(list1)+","+str(list2)+",1.17\")"
    os.chdir(Path)
    textfile=["global color=yellow dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n fk5"]+list(map(converter,radio_data[...,0],radio_data[...,1]))
    np.savetxt('radioRegionsFiltered'+str(obsID)+'.reg',textfile,newline='\n',fmt="%s")

#===================================================
#updating the acis evt file
def TCRVL(filename,offset,Path):
    #loads the fits hdulist 
    hdulist=fits.open(filename, mode='update')
    
     # Parse the WCS keywords in the primary HDU
    w = wcs.WCS(hdulist[0].header)
      
    TCRVL11=hdulist[1].header['TCRVL11'] 
    TCRVL12=hdulist[1].header['TCRVL12']
    
    if not os.path.exists(Path+'coordinate_flag.txt'):
        hdulist[1].header['TCRVL11']=hdulist[1].header['TCRVL11']+offset[0]
        hdulist[1].header['TCRVL12']=hdulist[1].header['TCRVL12']+offset[1]
    
        hdulist.flush()
        hdulist.close()

        print("Applying shift to event file")

    else:
        print("Coordinate change to event file  has already been done")

#======================================================
#updating the asol1 fits file
def asol1_update(asol1File,offset,Path):
    hdulist=fits.open(asol1File,mode='update')
    w=wcs.WCS(hdulist[0].header)

    RA_NOM=hdulist[1].header['RA_NOM'] 
    DEC_NOM=hdulist[1].header['DEC_NOM']
    
    if not os.path.exists(Path+'coordinate_flag.txt'):
        hdulist[1].header['RA_NOM']=hdulist[1].header['RA_NOM']+offset[0]
        hdulist[1].header['DEC_NOM']=hdulist[1].header['DEC_NOM']+offset[1]
    
        hdulist.flush()
        hdulist.close()

        print("Applying shift to asol1 file")
        f= open(Path+"coordinate_flag.txt","w+")
    else:
        print("Coordinate change to asol1 file has already been done")

#=======================================================
#converting sexa coordinates to degrees
def sexa_to_deg(sexa_coords):

#    from ChPipe_parfile.py import sgrApos,magpos
    flag=sexa_coords
    def isfloat(value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    array=[[],[],[]]
    y=0
    for i in sexa_coords:
        if i=='-' or isfloat(i) or i=='.':
            array[y]=np.append(array[y],i)
        else:
            y=y+1
    
    if flag==sgrApos[0] or flag==magpos[0]:
        deg_value=(float(''.join(array[0]))+ float(''.join(array[1]))/60+float(''.join(array[2]))/3600)*360/24

    elif flag==sgrApos[1] or flag==magpos[1]:
        deg_value=(float(''.join(array[0]))- float(''.join(array[1]))/60-float(''.join(array[2]))/3600)
        
    return deg_value
        
#========================================================
#copying the files

def copy_files(evt_file,asol_file,dataPath,obsID):

    
    wcsAsol = dataPath + str(obsID) + '_repro_wcs_asol1.fits'
    shutil.copy(asol_file,wcsAsol)
    wcsEvt = dataPath + str(obsID) + '_repro_wcs_evt2.fits'
    shutil.copy(evt_file, wcsEvt)

    return wcsEvt,wcsAsol
#+++++++++++++++++++++++++++++++++++++++
# ====================================== 
def wcsRoutines(obsID, squareEvtFile, fullevtFile, psList, wavsrcFile, dataPath, asol1File):

    print( "      starting wcs_match" )
    inFile = wavsrcFile
    refsrcFile = psList
    wcsFile = squareEvtFile
    xfmFile = dataPath + str(obsID) + 'xfm.fits'
    radius = 1
    residlim = 0.5
    verbose = 1
    clobber = 'yes'
    logfile = dataPath+str(obsID)+"_wcs_match.notes"
    rt.wcs_match(infile=inFile, refsrcfile=refsrcFile, wcsfile=wcsFile,
                 outfile=xfmFile, radius=radius, residlim=residlim,
                 verbose=verbose, clobber=clobber, logfile=logfile)
    
    print( "      starting wcs_update for asol file")
    transformFile = xfmFile
    wcsAsol = dataPath + str(obsID) + '_repro_wcs_asol1.fits'
    clobber = 'yes'
    rt.wcs_update.punlearn()
    rt.wcs_update(infile=asol1File, transformfile=transformFile, outfile=wcsAsol, clobber=clobber)

    print( "      starting wcs_update for evt file")
    wcsEvt = dataPath + str(obsID) + '_repro_wcs_evt2.fits'
    shutil.copy(fullevtFile, wcsEvt)
    transformFile = xfmFile
    rt.wcs_update.punlearn()
    rt.wcs_update(infile=wcsEvt, transformfile=transformFile, outfile='')

    return wcsEvt, wcsAsol

# ======================================
def makeLightCurves(srcreg, bkgreg, radtag, nametag, evtFile, obsID, dataPath, emin, emax, etag, lcbin, ccd):

    # define input file names
    srcFile = evtFile+'[sky=region('+srcreg+'),ccd_id='+str(ccd)+',energy='+str(emin)+':'+str(emax)+'][bin time=::'+str(lcbin)+']'
    srcFile2 = evtFile+'[sky=region('+srcreg+'),ccd_id='+str(ccd)+',energy='+str(emin)+':'+str(emax)+'][bin time=::0.44104]'
    bkgFile = evtFile+'[sky=region('+bkgreg+'),ccd_id='+str(ccd)+',energy='+str(emin)+':'+str(emax)+']'
    # define output file names
    lc_nobkg = str(obsID)+"_"+nametag+"_"+etag+"_"+radtag+"_lc_"+str(lcbin)+"s.fits"
    lc_nobkg_noTbin = str(obsID)+"_"+nametag+"_"+etag+"_"+radtag+"_lc_0.44104s.fits"
    lc_bkg   = str(obsID)+"_"+nametag+"_bkgsub_"+etag+"_"+radtag+"_lc_"+str(lcbin)+"s.fits"
    
    # source without background
    rt.dmextract(srcFile, outfile=dataPath+lc_nobkg, bkg="", opt="ltc1", clobber="yes")
    print( "      created "+lc_nobkg)
    rt.dmextract(srcFile2, outfile=dataPath+lc_nobkg_noTbin, bkg="", opt="ltc1", clobber="yes")
    print( "      created "+lc_nobkg_noTbin)
    # source with background
    rt.dmextract(srcFile, outfile=dataPath+lc_bkg, bkg=bkgFile, opt="ltc1", clobber="yes")
    print( "      created "+lc_bkg)

    return dataPath+lc_nobkg, dataPath+lc_bkg

# ======================================
def plotLCcurve(lcfits_nobkg, lcfits_bkg, imgname, dataPath):

    # for no background
    hdulist = fits.open(lcfits_nobkg)
    lc1 = hdulist[1].data
    photon_time1 = lc1.field('time')
    photon_counts1 = lc1.field('counts')
    # for with background
    hdulist = fits.open(lcfits_bkg)
    lc2 = hdulist[1].data
    photon_time2 = lc2.field('time')
    photon_counts2 = lc2.field('counts')

    fullimgname = dataPath+imgname+".jpg"
    plt.figure()
    plt.xlim( min(photon_time1), max(photon_counts1) )
    plt.plot(lc1.time, lc1.counts)
    plt.plot(lc2.time, lc2.counts)
    plt.ylabel("Count Rate")
    plt.xlabel("Time")
    plt.title(fullimgname)
    print( "      created lightcurve image "+ imgname)
    return

# ======================================
def plotCiaoLcurve(lcfits_nobkg, lcfits_bkg, nametag, lcbin, obsID, dataPath):    

    imgname = str(obsID)+"_"+nametag+"_lcurve_"+str(lcbin)+"s.png"

    # http://cxc.harvard.edu/chips/faq/script.html
    lcstr_nobkg = lcfits_nobkg+'[cols time_bin,count_rate]'
    lcstr_bkg = lcfits_bkg+'[cols time_bin,count_rate]'
    
    pc.set_preference("window.display", "false")
    hlui.make_figure(lcstr_nobkg, ["symbol.style","none"])
    hlui.add_curve(lcstr_bkg, ["symbol.style","none", "line.color","red"] )

    print( "      created lightcurve image ")
    print( "        --> "+dataPath+imgname)
    pc.print_window(dataPath+imgname, ["clobber", 'True'])
    

    # http://cxc.harvard.edu/chips/ahelp/chips4.html#Running_ChIPS_in_batch_mode
    #ChIPS always requires the presence of a X server, even when used in a script or
    #program which does not produce on-screen plots (e.g. the "window.display" preference
    #has been set to "false"). When a ChIPS server is started, it tries to connect to the
    #X-display indicated by the DISPLAY environment variable, which can result in one of three cases:
    #  DISPLAY environment variable is valid
    #    In this case the ChIPS server will use the indicated X display for its plots.
    #  DISPLAY environment variable is invalid or unset
    #    it will try to start up a temporary, virtual, framebuffer using the Xvfb program.
    #    If this works then the ChIPS commands will work, including calls to print_window(),
    #    but there will be no on-screen display.
    #  No DISPLAY and no Xvfb
    #    If the Xvfb program does not exist on your system then the server will fail to start
    #    up and report the following:
    #    chips ERROR: Failed to open DISPLAY or Xvfb session- exiting chipsServer

    return

# ====================================== 
def getSpectrum(evt2, srcReg, nametag, rtag, hkfiles, obsID, dataPath, bkgReg=''):

    rt.specextract.punlearn()
    [asol1, asol1lis, msk1, bpix] = hkfiles

    inEvtFile = evt2+"[sky=region("+srcReg+")]"
    root = dataPath+str(obsID)+"_"+nametag+"_"+rtag

    if bkgReg=='':
        bkg = ''
    else:
        bkg = evt2+"[sky=region("+bkgReg+")]"
    
    rt.specextract.punlearn()
    #print( inEvtFile)
    #print( asol1)
    #print( root)
    rt.specextract(infile=inEvtFile, outroot=root,  bkgfile=bkg, asp=asol1, mskfile=msk1, badpixfile=bpix, 
                   weight='no', correctpsf='yes', clobber='yes')
                   

# ======================================
if __name__ == "__main__":                 
    Run(sys.argv[1:])
