source ~psutton/opt/xpipeline/r5650_matlab2015a_localMCR/bin/envsetup.sh 
source activate /home/scott.coughlin/.conda/envs/xpipeline-py36/
xpipeline-workflow --search-type grb --params-file grb_mini.ini -n GRB160830 -g 1128477617 -r 307.65 -d 45.72  -i H1 L1 --FAP 99.99 --catalogdir /home/scoughlin/PhD/src/xpipeline-waveforms/  --off-source-inj
cp /home/scott.coughlin/PhD/src/xpipeline-waveforms/* input/
rm -rf output/off_source/ output/on_source/ output/ul_source/ output/simulations_*
