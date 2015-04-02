# make the tmplist of files to run over
sed $1'q;d' listfile.txt | sed 's/eoscms\/\/eos\/cms/xrootd-cms.infn.it\//' > tmplist.txt

# run RazorRun
./RazorRun tmplist.txt Analyses ntupleName.root 0

# run cmsRun to produce a valid FrameworkJobReport.xml 
cmsRun -j FrameworkJobReport.xml -p PSet.py
