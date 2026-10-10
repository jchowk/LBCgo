# LBCgo Gain/ReadNoise: 2014.06 calibration

# Environment for LBcgo
```
conda activate lbcgo
```

# Data directory: 
```
export LBCGO_RAW=/Users/howk/Dropbox/Data/LBT/Raw/2014.06/LBC/
```

# Create inputs.ecsv:
```
python ~/python/LBCgo/calibration/make_inputs.py $LBCGO_RAW -o inputs_full.ecsv --overwrite --levels
```

The final `inputs.ecsv` is an edited version of the original output removing saturated flats and extraneous biases. 

# Examine files in the directory
```
from ccdproc import ImageFileCollection
ic0b = ImageFileCollection('./', keywords=keywds,
                                 glob_include='lbcb*fits.gz')
ic0r = ImageFileCollection('./', keywords=keywds,
                                 glob_include='lbcr*fits.gz')
```

# Run run.py code:
```
LBCGO_RAW=/Users/howk/Dropbox/Data/LBT/Raw/2016.02/ python run.py --box 0    
```

This could just be 
```
export LBCGO_RAW=/Users/howk/Dropbox/Data/LBT/Raw/2016.02/
python run.py --box 0    
```


# Claude prompt to do some of this

```
Full init not needed. This directory is meant to run the gain and read noise calibrations for LBCgo using the 02/2016 dataset I have available (see the file `whatidid.md`). All work should be done in the Conda environment `lbcgo`. This is in support of the LBCgo development (https://github.com/jchowk/LBCgo ; see also the directories above this, notably /Users/howk/python/LBCgo/docs/planning/). 

I've run the `make_inputs.py` script for this dataset following through line 15 in the `whatidid.md` file. Can you examine `inputs_full.ecsv` to select calibration "sets" for deriving the gain and readnoise, putting the culled version in `inputs.ecsv`. 
```