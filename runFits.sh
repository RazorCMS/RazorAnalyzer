
for cat in DiJet_0b #DiJet_1b DiJet_2b #MultiJet_0b MultiJet_1b MultiJet_2b LeptonJet_0b LeptonJet_1b LeptonJet_2b LeptonMultiJet_0b LeptonMultiJet_1b LeptonMultiJet_2b #DiJet MultiJet LeptonJet LeptonMultiJet
do
#echo $cat
python python/RazorFitMacro.py --config config/run2_2017_03_13_SeparateBtagFits_forToys.config $cat #> ${cat}.txt
#python python/RazorFitMacro.py --config config/run2_2017_03_13_SeparateBtagFits_forToys.config --load-fit --run-toys $cat #> ${cat}.txt
done