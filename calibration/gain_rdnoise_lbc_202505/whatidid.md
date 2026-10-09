# Data in the directory: 

```
export LBCGO_RAW='/Users/howk/Dropbox/Data/LBT/Raw/2025.05_calib/'
```

# Create inputs.ecsv:
```
python ~/python/LBCgo/calibration/make_inputs.py ./ -o inputs_full.ecsv --overwrite --levels
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
LBCGO_RAW=/Users/howk/Dropbox/Data/LBT/Raw/2025.05_calib python run.py --box 0    
```

This could just be 
```
export LBCGO_RAW=/Users/howk/Dropbox/Data/LBT/Raw/2025.05_calib 
python run.py --box 0    
```