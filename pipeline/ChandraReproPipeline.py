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

# ==========================================
def Run(argv):

    print "================================================ "
    print "    Beginning Chandra reprocessing pipeline"
    print "================================================ "
    
    # --------- import parameters from parfile
    if type(argv) is str:
        parfile = argv
    if type(argv) is list:
        parfile = argv[0]
    
    parfilename = str(parfile.split("/")[-1:][0])
    par = imp.load_source(parfilename, parfile)
        
    print "The following parameter file was successfully uploaded:"
    print "     "+parfile
    parholder = {key: value for key, value in par.__dict__.iteritems() if not (key.startswith('__') or key.startswith('_'))}
    for key, value in parholder.iteritems():
        print str(key)+" = "+str(value)
    
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
        print "-----------------"
        print "beginning download process"
        print "-----------------"
        downloadData(par.obsID_list, par.downloadPath)
                             
    # creates the repro directory as repro_Nov17
    tstamp= datetime.datetime.now()
    reproDir = "repro_"+tstamp.strftime("%b")+tstamp.strftime("%y")
    
    for obsID in par.obsID_list:
        print "-----------------"
        print "beginning analysis for "+str(obsID)
        print "-----------------"
        
        # first process the requred full paths
        reproPath = par.reproPathBase + "/" + reproDir + "/" + str(obsID) + "/"
        dlPathFull = par.downloadPath+"/"+str(obsID)
        dlPrimPath = par.downloadPath+"/"+str(obsID)+"/primary/"

        print "** creating directory for reprocessing data"
        if not os.path.exists(reproPath):
                os.makedirs(reproPath)
        print "   "+reproPath

        if par.reproBool=='yes' or par.reproBool=='YES' or par.reproBool=="Yes":
            print "** reprocessing data"
            rt.chandra_repro(indir=dlPathFull, outdir=reproPath, root='acisf'+str(obsID), badpixel='yes',
                             process_events='yes', set_ardlib='yes', pix_adj='edser', check_vf_pha='no',
                             cleanup='yes', destreak='yes', clobbe='yes')
            print "   finished reprocessing"

        print "** finding the event, asol, mask, and badpixel files in reprecessed dataset"
        evt2File, asol1,  asol1lis, msk1, bpix  = findFiles(reproPath)   #contain the full paths and file names
        hkfiles = [asol1,  asol1lis, msk1, bpix]

        print "** setting the badpix file "
        setbpixFile(reproPath, bpix)

        print "** converting SrgrA, magnetar, and background  ra/dec into x/y coordinates"
        [sgrApix, magpix, backpix] = getPixelCoords(reproPath, evt2File, [par.sgrApos, par.magpos, par.bkgpos], asol1)
        
        print "** creating the energy- and spatial-cut event file for wavdetect"
        squareRegFile = makeRegFiles(reproPath, str(obsID)+"_wvdet_square.reg", sgrApix, [boxwidth, boxwidth], "box")
        squareEvtFile = dmcopyEvtFile(evt2File, reproPath, squareRegFile, obsID, "_square_2-8keV", erange=[emin, emax], ftype='evt')
        
        print '** starting getPSFmap'
        psfmap = getPSFmap(squareEvtFile)

        print "** starting wavdetect call"
        wavdetsrcFile = wvdetectCall(squareEvtFile, psfmap)
        
        print "** creating regions and ds9 image before wcs corrections"
        sgrAregFile = makeRegFiles(reproPath, str(obsID)+"_sgrA_orig.reg", sgrApix, par.sgrArad, "circ")
        magregFile = makeRegFiles(reproPath, str(obsID)+"_magnetar_orig.reg", magpix, par.magrad, "circ")
        getds9image([sgrAregFile, magregFile], reproPath+ str(obsID)+"_repro_2-8keV_prewcs.jpg", squareEvtFile)
        
        print "** performing wcs corrections"
        evt2wcsFile, asolwcs = wcsRoutines(obsID, squareEvtFile, evt2File, psListFile, wavdetsrcFile, reproPath, asol1)
        print "   "+evt2wcsFile
        hkfiles[0] = asolwcs

        print "** creating regions and ds9 image after wcs corrections"
        [sgrApixwcs, magpixwcs, backpixwcs] = getPixelCoords(reproPath, evt2wcsFile, [par.sgrApos, par.magpos, par.bkgpos], asol1)
        sqReg_wcsFile = makeRegFiles(reproPath, str(obsID)+"_wvdet_square_wcs.reg", sgrApixwcs, [boxwidth, boxwidth], "box")
        sgrAreg_wcsFile = makeRegFiles(reproPath, str(obsID)+"_sgrA_wcs.reg", sgrApixwcs, par.sgrArad, "circ")
        magreg_wcsFile = makeRegFiles(reproPath, str(obsID)+"_magnetar_wcs.reg", magpixwcs, par.magrad, "circ")
        backreg_wcsFile = makeRegFiles(reproPath, str(obsID)+"_bkg_wcs.reg", magpixwcs, par.bkgrad, "ann")
        cutEvtwcsFile = dmcopyEvtFile(evt2wcsFile, reproPath, sqReg_wcsFile , obsID, "_square_2-8keV_wcs", erange=[emin, emax], ftype='evt')
        getds9image([sgrAreg_wcsFile, magreg_wcsFile,backreg_wcsFile], reproPath+str(obsID)+"_repro_2-8keV_postwcs.jpg", cutEvtwcsFile)


        # -----------
        print "** beginning event file region extractions"
        sgrEcutwcsEvt = dmcopyEvtFile(evt2wcsFile, reproPath, sgrAreg_wcsFile , obsID, "_SgrA_2-8keV_wcs_evt", erange=[emin,emax], ftype='evt')
        sgrwcsEvt     = dmcopyEvtFile(evt2wcsFile, reproPath, sgrAreg_wcsFile , obsID, "_SgrA_wcs_evt", ftype='evt')
        magEcutwcsEvt = dmcopyEvtFile(evt2wcsFile, reproPath, magreg_wcsFile , obsID, "_magnetar_2-8keV_wcs_evt", erange=[emin,emax], ftype='evt')
        magwcsEvt     = dmcopyEvtFile(evt2wcsFile, reproPath, magreg_wcsFile , obsID, "_magnetar_wcs_evt", ftype='evt')
        
        print "** beginning the light curve extractions for 2-8 keV"
        # the outputs here are lc file names"
        sgr_LCnobkg, sgr_LCbkg = makeLightCurves(sgrAreg_wcsFile, backreg_wcsFile, sradtag, "sgra", sgrEcutwcsEvt, obsID, reproPath,
                                                 emin, emax, etag, par.lcbin, par.ccd_id)
        mag_LCnobkg, mag_LCbkg = makeLightCurves(magreg_wcsFile, backreg_wcsFile,  mradtag, "magnetar", magEcutwcsEvt, obsID, reproPath,
                                                 emin, emax, etag, par.lcbin, par.ccd_id)
        print "** plotting lightcurves"
        plotCiaoLcurve(sgr_LCnobkg, sgr_LCbkg, "sgrA", par.lcbin, obsID, reproPath)
        plotCiaoLcurve(mag_LCnobkg, mag_LCbkg, "magnetar", par.lcbin, obsID, reproPath)
        
        
        print "** creating spectral file"        
        sgrSpecFile = getSpectrum(evt2wcsFile, sgrAreg_wcsFile, "sgra", sradtag, hkfiles, obsID, reproPath, bkgReg=backreg_wcsFile)
        magSpecFile = getSpectrum(evt2wcsFile, magreg_wcsFile, "magnetar", mradtag, hkfiles, obsID, reproPath, bkgReg=backreg_wcsFile)
        
        
# ======================================
# ======================================
def downloadData(idList, downloadPath):
        
    currentDir = os.getcwd()
    os.chdir(downloadPath)
    print "downloading data to the download dir"
    print "    --> "+downloadPath

    # this controls verbosity of download
    lw.initialize_logger('download', verbose=1)
    # Downloads this to the download (local) directory
    for idNumb in idList:
        download_chandra_obsids([idNumb])
        # if needed, unzip the files
        print "unzipping"
        os.system("gunzip -r "+downloadPath+"/"+str(idNumb) )
        
    # change back to the original directory    
    os.chdir(currentDir)

    
# ======================================
def findFiles(dataPath):

    foundfiles = glob.glob(dataPath + '*evt2.fits')
    if len(foundfiles) == 1:
        evt2file = foundfiles[0]
    else:
        print 'error: either none or many evt2 files have been found'
        quit()
    foundfiles = glob.glob(dataPath + '*asol1.fits')
    if len(foundfiles) == 1:
        asol1file = foundfiles[0]
    else:
        print 'error: either none or many asol1.fits files have been found'
        quit()
    foundfiles = glob.glob(dataPath + '*asol1.lis')
    if len(foundfiles) == 1:
        asol1lis = foundfiles[0]
    else:
        print 'error: either none or many asol1.lis files have been found'
        quit()
    foundfiles = glob.glob(dataPath + '*msk1.fits')
    if len(foundfiles) == 1:
        msk1file = foundfiles[0]
    else:
        print 'error: either none or many msk files have been found'
        quit()
    foundfiles = glob.glob(dataPath + '*repro_bpix1.fits')
    if len(foundfiles) == 1:
        bpixfile = foundfiles[0]
    else:
        print 'error: either none or many repro_bpix1 files have been found'
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
        print "   pix     "+str(pix[0])+", "+str(pix[1])
        outArray.append(pix)

    return outArray
    '''
    rt.dmcoords(evt2, asolfile=asol1, ra=sPos[0], dec=sPos[1], option='cel', verbose=1)
    sPix = [rt.dmcoords.x, rt.dmcoords.y]
    print "   sgrA pix     "+str(sPix[0])+", "+str(sPix[1])
    rt.dmcoords(evt2, asolfile=asol1, ra=mPos[0], dec=mPos[1], option='cel', verbose=1)
    mPix = [rt.dmcoords.x, rt.dmcoords.y]
    print "   magnetar pix "+str(mPix[0])+", "+str(mPix[1])
    rt.dmcoords(evt2, asolfile=asol1, ra=bPos[0], dec=bPos[1], option='cel', verbose=1)
    bPix = [rt.dmcoords.x, rt.dmcoords.y]
    print "   bkg pix      "+str(bPix[0])+", "+str(bPix[1])
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
    print "      created "+regFile
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
        
    print "      created "+outEvtName
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

# ====================================== 
def wcsRoutines(obsID, squareEvtFile, fullevtFile, psList, wavsrcFile, dataPath, asol1File):

    print "      starting wcs_match" 
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
    
    print "      starting wcs_update for asol file"
    transformFile = xfmFile
    wcsAsol = dataPath + str(obsID) + '_repro_wcs_asol1.fits'
    clobber = 'yes'
    rt.wcs_update.punlearn()
    rt.wcs_update(infile=asol1File, transformfile=transformFile, outfile=wcsAsol, clobber=clobber)

    print "      starting wcs_update for evt file"
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
    print "      created "+lc_nobkg
    rt.dmextract(srcFile2, outfile=dataPath+lc_nobkg_noTbin, bkg="", opt="ltc1", clobber="yes")
    print "      created "+lc_nobkg_noTbin
    # source with background
    rt.dmextract(srcFile, outfile=dataPath+lc_bkg, bkg=bkgFile, opt="ltc1", clobber="yes")
    print "      created "+lc_bkg

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
    print "      created lightcurve image "+ imgname
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

    print "      created lightcurve image "
    print "        --> "+dataPath+imgname
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
    #print inEvtFile
    #print asol1
    #print root
    rt.specextract(infile=inEvtFile, outroot=root,  bkgfile=bkg, asp=asol1, mskfile=msk1, badpixfile=bpix, 
                   weight='no', correctpsf='yes', clobber='yes')
                   

# ======================================
if __name__ == "__main__":                 
    Run(sys.argv[1:])
