#manual trans-rot script

import numpy as np
from astropy.io import fits
from astropy import wcs
import sys
import csv

import ciao_contrib.runtool as rt
import glob
import imp
from subprocess import call
import os
import datetime
import ciao_contrib.logger_wrapper as lw
from ciao_contrib.cda.data import download_chandra_obsids
import shutil
import pychips as pc
import pychips.hlui as hlui
#==============================================================================================

#IMPORTANT PARAMETERS TO FILL IN BEFORE EXECUTING SCRIPT

#==============================================================================================

#radio coordinates of sag A* and magnetar

mag_rad_ra=266.4173708
mag_rad_dec=-29.008289

sag_rad_ra=266.416837
sag_rad_dec=-29.0078105

obsID=15043
downloadPath='/data/irulan/gc_catalog/students/alexk/chandra_obsID_list'
Script_Path='/data/irulan/gc_catalog/students/alexk/python_scripts'
RegFilePath='/data/irulan/gc_catalog/students/alexk/regionfiles'
PATH =downloadPath+ '/'+str(obsID)+'/repro'
tool=False
all_sources_average_translation=True

#assigning the files to certain names
file_name= 'cdfs'+str(obsID)+'_broad_thresh.img'
acisf_file='acisf'+str(obsID)+'_repro_evt2.fits'


#==============================================================================================


# ALL THE FUNCTIONS ARE LISTED BELOW THIS POINT


#==============================================================================================



def deg_to_arc(x):
    y=x*3600
    return y 

def getFluximage(obsID):
    infile= 'acisf'+str(obsID)+'_repro_evt2.fits'
    outroot= 'cdfs'+str(obsID)
    bands='broad'
    bin=1
    
    rt.fluximage.punlearn()
    rt.fluximage(infile=infile, outroot=outroot, bin=bin, bands=bands)

def getPSFmap(filename):
    
    psfmap = 'cdfs'+str(obsID) + '_psfmap.fits'
    rt.mkpsfmap(infile=filename, outfile=psfmap, energy=3.8, ecf=0.393, clobber='yes')
    return psfmap

def CRVAL(filename):

    #loads the fits hdulist 
    hdulist=fits.open(filename, mode='update')
    
     # Parse the WCS keywords in the primary HDU
    w = wcs.WCS(hdulist[0].header)

    # Print out the RA and DEC position of the tangent point
    
    CRVAL1= w.wcs.crval[0]
    CRVAL2 =w.wcs.crval[1]

    return CRVAL1,CRVAL2

#Running wavdetect
def getWavdetect(filename):

    infile= filename
    outFile = filename[:-11] + '_src.fits'
    scellFile = filename[:-11] + '_scell.fits'
    imageFile = filename[:-11] + '_img.fits'
    nbkgFile = filename[:-11] + '_nbkg.fits'
    regFile = filename[:-11] + '_src.reg'
    scales = '1.0 2.0 4.0 8.0 16.0'
    sigThresh = 1e-06
    psfmap= 'cdfs'+str(obsID) + '_psfmap.fits'
    expfile= file_name[:-3]+'expmap'
    rt.wavdetect.punlearn()
    rt.wavdetect(infile=infile, psffile=psfmap, outfile=outFile,
                 scellfile=scellFile, imagefile=imageFile,
                 defnbkgfile=nbkgFile, regfile=regFile,
                 scales=scales, sigthresh=sigThresh, expfile=expfile, clobber='yes')

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

def rotation(phi):

    
    total_angle= -(np.pi/2)+theta3+phi
    
    return total_angle

def scaling():
    
    print("SCALING")
    distance_rad=((mag_rad_ra-sag_rad_ra)**2+(mag_rad_dec-sag_rad_dec)**2)**(0.5)
    scalefac=distance_rad/L3
    return scalefac

def wcsRoutines(filename):
    
    infile = filename
    wcsEvt = 'acisf'+str(obsID)+'_repro_evt2.fits'
    rotang= rotation_angle+0.8
    #the number I am adding isn't legit, the goal_phi and the detected phi have to match for it to be legit
    outfile= filename[:-11]+'2_thresh.img'
    clobber='yes'
    
    print("Copying thresh file")
    rt.dmcopy(infile=infile,outfile=outfile ,clobber=clobber)

    print("starting wcs_update for evt file")
    rt.wcs_update.punlearn()
    rt.wcs_update(infile=outfile ,wcsfile=wcsEvt, rotang=rotang,scalefac=scalefac, outfile='')
    print("Done")


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
        

    while not while_var and counts+1<len(allsources):
        counts=counts+1
        while_var=in_circle(sag_rad_ra,sag_rad_dec,allsources[counts][0],allsources[counts][1],2/3600)[0]
        

    sag_x_ra=in_circle(sag_rad_ra,sag_rad_dec,allsources[counts][0],allsources[counts][1],2/3600)[1]
    sag_x_dec=in_circle(sag_rad_ra,sag_rad_dec,allsources[counts][0],allsources[counts][1],2/3600)[2]
    print("Sag A* wavdetect output")
    print(sag_x_ra,sag_x_dec)

    
    ra_mag_offset=mag_rad_ra-src[0][0]
    dec_mag_offset=mag_rad_dec-src[0][1]

    ra_sag_offset=sag_rad_ra-sag_x_ra
    dec_sag_offset=sag_rad_dec-sag_x_dec

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
 
def getDmcopy(filename,outfile):

    infile=filename
    outfile=outfile
    rt.dmcopy(infile=infile,outfile=outfile ,clobber='yes')


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

"""
==============================================================================================================================================================
Loading regions from ascii radio file and filtetering process beginning


==============================================================================================================================================================
"""

def load_regions():
    os.chdir(Script_Path)


    with open('chandra_point_sources_wsagmag.ascii', 'rt') as f:
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
    #---------------------------------------------
    #comparing coordinates to see if they match

"""
===================================================================
This function takes the outputs from the wavdetect source file and
filters the radio region catalogue file in the input. 

===================================================================
"""


def filter_regions():

    
    def in_circle(center_x, center_y, x, y,radius):

        global filtered_rad_reg
        global wav_regions
        square_dist = (center_x - x) ** 2 + (center_y - y) ** 2
        if square_dist <= radius ** 2:

            filtered_rad_reg=np.concatenate((filtered_rad_reg,np.array([[x,y]])))

            wav_regions=np.array(np.append(wav_regions,i))
          
        return filtered_rad_reg

    for i in range (0,len(all_sources[...,0])):
            
       for j in range (0,len(reg[:][0])):
            coords=[all_sources[i][0],all_sources[i][1],reg[0][j],reg[1][j],1/3600]
            filtered_rad_reg=in_circle(*coords)
            
    filtered_rad_reg=np.delete(filtered_rad_reg,0,0)        
    return filtered_rad_reg

def all_translation():
    mean_offset=np.array([np.mean(filtered_radio_data[...,0]-filtered_wav_data[...,0]),np.mean(filtered_radio_data[...,1]-filtered_wav_data[...,1])])

    #if i install pandas, then the following command will be executable to make nice tables
    #pd.DataFrame([(deg_to_arc(filtered_radio_data[...,0]-filtered_wav_data[...,0])),(deg_to_arc(filtered_radio_data[...,1]-filtered_wav_data[...,1]))],["RA offset","DEC offset"])

    print("RA offset for each source")
    print(deg_to_arc(filtered_radio_data[...,0]-filtered_wav_data[...,0]))
    print("DEC offset for each source")
    print(deg_to_arc(filtered_radio_data[...,1]-filtered_wav_data[...,1]))
    print("The average offset is: ")
    print(deg_to_arc(mean_offset))
    return mean_offset


def save_reg():
    global filtered_radio_data
    def converter(list1,list2):
        return "circle("+str(list1)+","+str(list2)+",1.17\")"
                                                                                                                                                                                       
    textfile=["global color=yellow dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n fk5"]+list(map(converter,filtered_radio_data[...,0],filtered_radio_data[...,1]))

    os.chdir(RegFilePath)

    np.savetxt('radioRegionsFiltered'+str(obsID)+'.reg',textfile,newline='\n',fmt="%s")

    os.chdir(PATH)
    np.savetxt('radioRegionsFiltered'+str(obsID)+'.reg',textfile,newline='\n',fmt="%s")

"""
==================================================================================================
Filtering process ended. All has been saved.


===================================================================================================
"""


#======================================================================
#updating the acis evt file
def TCRVL(filename):
    #loads the fits hdulist 
    hdulist=fits.open(filename, mode='update')
    
     # Parse the WCS keywords in the primary HDU
    w = wcs.WCS(hdulist[0].header)
      
    TCRVL11=hdulist[1].header['TCRVL11'] 
    TCRVL12=hdulist[1].header['TCRVL12']
    
    if not os.path.exists('coordinate_flag.txt'):
        hdulist[1].header['TCRVL11']=hdulist[1].header['TCRVL11']+mean_offset[0]
        hdulist[1].header['TCRVL12']=hdulist[1].header['TCRVL12']+mean_offset[1]
    
        hdulist.flush()
        hdulist.close()

        print("Applying shift to event file")

    else:
        print("Coordinate change to event file  has already been done")

def asol1_update(asol1File):
    hdulist=fits.open(asol1File,mode='update')
    w=wcs.WCS(hdulist[0].header)

    RA_NOM=hdulist[1].header['RA_NOM'] 
    DEC_NOM=hdulist[1].header['DEC_NOM']
    
    if not os.path.exists('coordinate_flag.txt'):
        hdulist[1].header['RA_NOM']=hdulist[1].header['RA_NOM']+mean_offset[0]
        hdulist[1].header['DEC_NOM']=hdulist[1].header['DEC_NOM']+mean_offset[1]
    
        hdulist.flush()
        hdulist.close()

        print("Applying shift to asol1 file")
        f= open("coordinate_flag.txt","w+")
    else:
        print("Coordinate change to asol1 file has already been done")
        
    

    
def getDs9(filename,filename2):


    # this creates a command to pass to the os.
    # creates a jpg image that overlays the regions on the evt File image

    ds9Routine = "ds9 "+filename+" "
    
    ds9Routine = ds9Routine + " -zoom 8 -scale log -cmap b"
    ds9Routine= ds9Routine + " -pan to "+ str(mag_rad_ra)+ " " + str(mag_rad_dec)+ " wcs" 
    #ds9Routine = ds9Routine + " -regions "+ filename[:-11] + "_src.reg"
    
    ds9Routine= ds9Routine+ " "+ filename2
    ds9Routine= ds9Routine + " -pan to "+ str(mag_rad_ra)+ " " + str(mag_rad_dec)+ " wcs"
    ds9Routine= ds9Routine+ " -regions load all "+ "radioRegionsFiltered"+str(obsID)+".reg"
    ds9Routine = ds9Routine + " &"
    # append the image name "
    # call ds9 to the command line
    os.system(ds9Routine)

#================================================================================================

#RUNNING ALL THE FUNCTIONS

#================================================================================================
if not os.path.isdir(PATH[:-6]):
    os.chdir(downloadPath)
    #downloading the obs_ID if it isn't already downloaded
    download_chandra_obsids([obsID])
else:
    print("This obsID has already been downloaded")
#reprocess data
if not os.path.isdir(PATH):

    print ("** reprocessing data")
    os.chdir(PATH[:-6])    
    reproPath=PATH
    rt.chandra_repro.punlearn()
    rt.chandra_repro(indir=PATH[:-6], outdir=reproPath, root='acisf'+str(obsID), badpixel='yes',process_events='yes', set_ardlib='yes', pix_adj='edser', check_vf_pha='no',cleanup='yes', destreak='yes', clobbe='yes')
    print ("   finished reprocessing")

#locate yourself in the directory in which you want the files to download

os.chdir(PATH)

#making the Fluximage
print("Making Fluximage")
if not os.path.exists(file_name):
    getFluximage(obsID)
    print("Done Fluximage")
else:
    print("Fluximage already exists")


#making the PSF map    
print("Making PSFmap")
if not os.path.exists('cdfs'+str(obsID) + '_psfmap.fits'):
    getPSFmap(file_name)
    print("Done PSFmap")
else:
    print("PSFmap already exists")

#extrcting the CRVAL values from the image files

CRVAL=CRVAL(file_name)

print("The RA coordinate of the tangent point is "+ str(CRVAL[0]))
print("The DEC coordinate of the tangent point is "+ str(CRVAL[1]))

#running wavdetect for initial image

print("Running wavdetect")

if not os.path.exists(file_name[:-11] + '_src.fits'):
    getWavdetect(file_name)
    print("Done wavdetect")
else:
    print("wavdetect files already exist")
    
#Finding the asol1 File

foundfiles = glob.glob(PATH + '/*asol1.fits')
if len(foundfiles) == 1:
    asol1file = foundfiles[0]
else:
    print( 'error: either none or many asol1.fits files have been found')
    quit()


#getting the source list from wavdetect
source_list=wavdetect_source(file_name[:-11]+'_src.fits')[0]
print(source_list)
all_sources=wavdetect_source(file_name[:-11]+'_src.fits')[1]
    

#------------------------------------------------------------------------


#if i use the rotation manual code, here is some of the calculations that are not well detailed and this part needs huge improvement

#calculating angle between sag and mag, and between sag and tangent point


#angle between mag and sagA* according to radio coordinates
goal_phi= np.arctan((mag_rad_ra-sag_rad_ra)/(mag_rad_dec-sag_rad_dec))

#angle between mag and sagA* according to wavdetect
wav_phi= np.arctan((source_list[0][0]-source_list[2][0])/(source_list[0][1]-source_list[2][1]))

L1=((CRVAL[0]-source_list[0][0])**2+(CRVAL[1]-source_list[0][1])**2)**(0.5)
L2=((CRVAL[0]-source_list[2][0])**2+(CRVAL[1]-source_list[2][1])**2)**(0.5)
L3=((source_list[0][0]-source_list[2][0])**2+(source_list[0][1]-source_list[2][1])**2)**(0.5)

theta3= np.arccos((-(L1**2)+L2**2+L3**2)/(2*L2*L3))


#end of need to improve part
#---------------------------------------------------------------------------

#calling the rotation --need to improve (not actually performing the rotation) 

rotation_angle=(rotation(goal_phi)-rotation(wav_phi))

#calling the scalling function to get the scale factor (integer)

scalefac=scaling()

#defining the second file that is gotten using dmcopy to perform modifications on it 
file_name2='cdfs'+str(obsID)+'_broad2_thresh.img'


###This is the script for rotation and translation also. If ever we want to use this part of the script, we need to adjust the rotation angle calculations.  

if tool:
    wcsRoutines(file_name)
    file_name2='cdfs'+str(obsID)+'_broad2_thresh.img'

    print("Running wavdetect for updated file")
    getWavdetect(file_name2)
    print("Done wavdetect")

    source_list2=wavdetect_source(file_name2[:-11]+'_src.fits')[0]
    all_sources2=wavdetect_source(file_name2[:-11]+'_src.fits')[1]
    wav_phi_new= np.arctan((source_list2[0][0]-source_list2[2][0])/(source_list2[0][1]-source_list2[2][1]))

    print("Goal angle: "+ str(goal_phi))
    print("Detected angle: "+ str(wav_phi))
    print("New angle: "+ str(wav_phi_new))
    print("Scaling factor: "+ str(scalefac))

#exectuing translation code

if  tool:
    trans=translation(source_list2,all_sources2)
else:
    trans=translation(source_list,all_sources)

#copying the file again in order to perform the translation now (so we don't mess with the original file and can have a before and after comparaison

if  tool:
    getDmcopy(file_name2,file_name2[:-12]+'3_thresh.img')
    
else:
    getDmcopy(file_name,file_name[:-11]+'2_thresh.img')
    

#defining the new filename after having used dmcopy

file_name3='cdfs'+str(obsID)+'_broad3_thresh.img'
    
#loading the regions from the ascii radio file and defining variables that are used in the function

reg=load_regions()
wav_regions=np.zeros((1,1))
filtered_rad_reg=np.zeros((1,2))

#filtering the radio region data file and assigning them to variables, filtering the wavdetect regions as well

filtered_radio_data=filter_regions()    
wav_regions=np.ndarray.tolist(np.delete(wav_regions,0,0))


filtered_wav_data=np.column_stack([np.take(all_sources[...,0],wav_regions),np.take(all_sources[...,1],wav_regions)])
print("On the left is the radio coordinates, and on the right are the wavedetect x-ray coordinates.")
print(np.column_stack([filtered_radio_data,filtered_wav_data]))

#Taking the offset of the distance of all sources

mean_offset=all_translation()

#saving the region files into the current directory and in the specified directory

save_reg()

#changing paths for this part since save_reg changes the directory to another specified directory

os.chdir(PATH)

#Applying the translation to the coordinate map using ciao dmhedit tool.
#assigning the values of the shift.
#if boolean all_sources_average_translation is set to true, the coordinate map is shifted by the average value of all source offset.


if  tool and not all_sources_average_translation:
    getDmhedit(file_name3,'CRVAL1',CRVAL[0]+trans[0])
    getDmhedit(file_name3,'CRVAL2',CRVAL[1]+trans[1])    
elif not tool and not all_sources_average_translation:
    getDmhedit(file_name2,'CRVAL1',CRVAL[0]+trans[0])
    getDmhedit(file_name2,'CRVAL2',CRVAL[1]+trans[1]) 
elif not tool and all_sources_average_translation:
    getDmhedit(file_name2,'CRVAL1',CRVAL[0]+mean_offset[0])
    getDmhedit(file_name2,'CRVAL2',CRVAL[1]+mean_offset[1]) 
    

if  tool:    
    getWavdetect(file_name3)
    source_list3=wavdetect_source(file_name3[:-11]+'_src.fits')[0]
    all_sources3=wavdetect_source(file_name3[:-11]+'_src.fits')[1]
else:
    getWavdetect(file_name2)
    source_list2=wavdetect_source(file_name2[:-11]+'_src.fits')[0]
    all_sources2=wavdetect_source(file_name2[:-11]+'_src.fits')[1]


#modifying the acis event file 

TCRVL=TCRVL(acisf_file)

#modifying the asol file

asol1_update(asol1file)


#having a before and after image using ds9

print("The new Transformation offset")
if tool:
    translation(source_list3,all_sources3)

else:
    translation(source_list2,all_sources2)

if tool:
    getDs9(file_name,file_name3)

else:
    getDs9(file_name,file_name2)

#other functions that I used but not in this script
"""
def sexa_to_deg(sexa_coords):

    def isfloat(value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    array=[[],[],[]]
    y=0
    for i in sexa_coords[0]:
        if isfloat(i) or i=='.':
            array[y]=np.append(array[y],i)
        else:
            y=y+1

    if sexa_coords==sgrApos[0] or sexa_coords==magpos[0]:
        deg_value=(float(''.join(array[0]))+ float(''.join(array[1]))/60+float(''.join(array[2]))/3600)*360/24

    elif sexa_coords==sgrApos[1] or sexa_coords==magpos[1]:
        deg_value=(float(''.join(array[0]))+ float(''.join(array[1]))/60+float(''.join(array[2]))/3600)
        
    return deg_value
"""
