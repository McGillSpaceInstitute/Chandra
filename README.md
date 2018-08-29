This repository contains a pipeline that does data analysis on sources in the galactic center. 
In this repository, we can find the pipeline folder.

In the pipeline folder, there are a certain set of documents: 
1.ChPipe_parfile.py
2.ChPipe_parfile.pyc
3.ChandraReproPipeline.py
4.chandra_point_sources_wmag.ascii
5.chandra_point_sources_wsagmag.ascii
6.manual_change.py

We can choose to completely ignore the .pyc file (#2).

The ChandraReproPipeline file is the main file which does all the analysis. In order to run the file,
you need to call ChPipe_parfile.py which call certain variables necessary for the execution of the file. 
Do not forget to change the variables in the parfile to match paths in your working directory. 

The chandra_point_sources_wmag.ascii file is a file that contains the RA and DEC coordinates of radio sources in the galactic center including the magnetar. 

The chandra_point_sources_wsagmag.ascii file is a file that contains the RA and DEC coordinates of radio sources in the galactic center including the magnetar and sagittarius A*. 

The manual_change.py file can be run on its own and matches X-ray coordinates to radio position coordinates in the galactic center. 
