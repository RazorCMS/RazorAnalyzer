### Put all relevant MADD plots into a pdf file for great justice

### Author: 
### Dustin Anderson

import os
import sys
import subprocess
from datetime import date

class Slide(object):
    def __init__(self, slide_title, files):
        """Arguments:
             slide_title: big title for slide
             files: list of file names"""
        self.slide_title = slide_title
        self.files = files
        self.width = 1.0/len(self.files)
        if self.width == 1.0:
            self.width = 0.8 # let's not get too crazy here

    def write(self, f):
        """f: output file object"""
        if not self.check_files: return
        f.write("\\begin{frame} \n")
        f.write("\\begin{center} \n")
        f.write("\\frametitle{%s} \n" % self.slide_title)
        for path in self.files:
            f.write("\\includegraphics[width=%f\\paperwidth]{%s} \n" % (self.width,path))
        f.write("\\end{center} \n")
        f.write("\\end{frame} \n")

    def append_dir_prefix(self, prefix):
        self.files = [ prefix+'/'+path for path in self.files ]

    def check_files(self):
        for f in self.files:
            if not os.path.isfile(f): 
                print "Didn't find file",f
                return False
        return True


def def_slides(plot_dir):
    """Returns: list of slide objects"""
    slides = []
    titles = { "TTJetsSingleLepton":"$t\\bar{t}$+jets control sample",
               "WJetsSingleLepton":"W+jets control sample",
               "WJetsSingleLeptonInv":"W+jets invisible control sample",
               }
    for control_region in ["TTJetsSingleLepton", "WJetsSingleLepton"]:
        slides += [
                Slide(titles[control_region], 
                    ["%s/MR_%s.pdf" % (control_region,control_region),
                     "%s/Rsq_%s.pdf" % (control_region,control_region)]),
                Slide(titles[control_region], 
                    ["%s/MRRsq_%sUnrolledDataMC.pdf" % (control_region,control_region),
                     "%s/%sScaleFactors.pdf" % (control_region,control_region.replace('SingleLepton',''))]),
                Slide(titles[control_region], 
                    ["%s/NBJetsMedium_%s.pdf" % (control_region,control_region),]),
                Slide(titles[control_region], 
                    ["%s/HT_%s.pdf" % (control_region,control_region),
                     "%s/MET_%s.pdf" % (control_region,control_region)]),
            ]
    control_region = "WJetsSingleLeptonInv"
    slides += [
            Slide(titles[control_region], 
                ["%s/MR_NoW_%s.pdf" % (control_region,control_region),
                 "%s/Rsq_NoW_%s.pdf" % (control_region,control_region)]),
            Slide(titles[control_region], 
                ["%s/MR_NoWRsq_NoW_%sUnrolledDataMC.pdf" % (control_region,control_region),
                 "%s/%sScaleFactors.pdf" % (control_region,control_region.replace('SingleLepton',''))]),
            Slide(titles[control_region], 
                ["%s/NBJetsMedium_%s.pdf" % (control_region,control_region),]),
        ]
    slides += [
            Slide("1L N$_{jets}$ correction", 
                ["OneLeptonForNJets/NJets40_SingleLepton.pdf"]),
            ]
    for jets in ['DiJet','MultiJet']:
        slides += [
                Slide("1L %s region, after corrections" % jets.lower(), 
                    ["OneLepton%sClosureTest/MRRsq_SingleLeptonUnrolledDataMC.pdf" % jets,
                     "OneLepton%sClosureTest/NBJetsMedium_SingleLepton.pdf" % jets,
                        ]),
                ]
        for B in ['0','1','2','3']:
            if jets == 'DiJet' and B == '3': continue
            slides += [
                Slide("1L %s region, after corrections -- %s B" % (jets.lower(),B),
                    ["OneLepton%sClosureTest%sB/MR_SingleLepton.pdf" % (jets,B),
                     "OneLepton%sClosureTest%sB/Rsq_SingleLepton.pdf" % (jets,B)])
                    ]
    photon_dir = "PhotonControlRegionMacro"
    photon_title = "photon+jets control region"
    slides += [
            Slide(photon_title,
                ["%s/Razor_PhotonControlRegion_PhotonPt_PhotonCR_Logy.pdf" % photon_dir,
                 "%s/Razor_PhotonControlRegion_PhotonEta_PhotonCR_Logy.pdf" % photon_dir]),
            Slide(photon_title,
                ["%s/Razor_PhotonControlRegion_MR_PhotonCR_Logy.pdf" % photon_dir,
                 "%s/Razor_PhotonControlRegion_Rsq_PhotonCR_Logy.pdf" % photon_dir]),
            Slide(photon_title,
                ["%s/Razor_PhotonControlRegion_MRRsqUnrolled_PhotonCR_Logy.pdf" % photon_dir,
                 "%s/Razor_PhotonControlRegion_NBtags_PhotonCR_Logy.pdf" % photon_dir]),
                ]
    tt2l_dir = "TTBarDileptonControlRegionMacro"
    for jets in ['DiJet','MultiJet']:
        tt2l_title = "t$\\bar{t}$ dilepton %s control region" % jets.lower()
        slides += [
                Slide(tt2l_title,
                    ["%s/Razor_TTBarDileptonCrossCheckRegion_MR_%s_Logy.pdf" % (tt2l_dir,jets),
                     "%s/Razor_TTBarDileptonCrossCheckRegion_Rsq_%s_Logy.pdf" % (tt2l_dir,jets)]),
                Slide(tt2l_title,
                    ["%s/Razor_TTBarDileptonCrossCheckRegion_MRRsqUnrolled_%s_Logy.pdf" % (tt2l_dir,jets),
                     "%s/Razor_TTBarDileptonCrossCheckRegion_NBtags_%s_Logy.pdf" % (tt2l_dir,jets)]),
                Slide(tt2l_title,
                    ["%s/Razor_TTBarDileptonCrossCheckRegion_NJets40_%s_Logy.pdf" % (tt2l_dir,jets),
                     "%s/Razor_TTBarDileptonCrossCheckRegion_MET_%s_Logy.pdf" % (tt2l_dir,jets)]),
                    ]


    for slide in slides:
        slide.append_dir_prefix(plot_dir)
    return slides

if __name__ == '__main__':
    if len(sys.argv) < 2:
        in_dir = "Plots/Razor2016"
    else:
        in_dir = sys.argv[1]

    fname = "MADD_plots_%s.tex"%date.today()
    with open(fname,'w') as f:
        f.write("\\documentclass{beamer}\n")
        f.write("\\beamertemplatenavigationsymbolsempty\n")
        f.write("\\usepackage{graphicx}\n")
        f.write("\\setbeamersize{text margin left=0pt,text margin right=0pt}\n")
        f.write("\\setbeamertemplate{footline}[frame number]\n")
        f.write("\\begin{document}\n")
        for slide in def_slides(in_dir):
            slide.write(f)
        f.write("\\end{document}\n")

    subprocess.call(['pdflatex', fname])
