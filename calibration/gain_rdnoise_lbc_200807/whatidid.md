# LBCgo Gain/ReadNoise: 2008.07 calibration

# Environment for LBCgo
```
conda activate lbcgo
```

# Data
Only the calibration frames downloaded from the LBT archive exist, as four tar files (observing log: `20080705.log.txt` in the same Raw directory):
```
~/Dropbox/Data/LBT/Raw/2008.07/Calibrations/IA2_LBT-SDT_12092008_*.tar
```
Extracted (outside Dropbox) with:
```
mkdir -p $LBCGO_RAW && cd $LBCGO_RAW
for f in ~/Dropbox/Data/LBT/Raw/2008.07/Calibrations/*.tar; do tar -xf "$f"; done
export LBCGO_RAW=<that directory>
```
What the tars contain (37 of the 96 calibration frames in the log): 10 LBCB and 10 LBCR biases (`25Bias_Bino`, 06:29-06:34 UT); 5 LBCB `SkyFlat_Usr_rot180` flats (SDT_Uspec, 1.48 s, 24-41k ADU); 5 LBCR i-SLOAN `SkyFlat_Vi_rot180_` flats (10.2 s, 22-40k ADU); 5 LBCR r-SLOAN `SkyFlat_Usr_rot180` flats (0.54 s, 42-65k ADU, saturated or above 0.7 x SATURATE); two LBCR test flats with a single extension (partial readouts, skipped by `make_inputs.py`). Missing: all LBCB V/B-BESSEL flats and the LBCR R-BESSEL flats.

# Create inputs
```
python ~/python/LBCgo/calibration/make_inputs.py $LBCGO_RAW -o inputs_full.ecsv --overwrite --levels
```
`inputs.ecsv`: four consecutive flat pairs per channel (LBCB SDT_Uspec; LBCR i-SLOAN), two biases each (first two biases per channel skipped). Adjacent pairs share a frame because only five usable flats exist per channel.

# Run
```
python run.py --box 0
```
(log of the run: scratchpad `run_2008.log`; not kept.)
