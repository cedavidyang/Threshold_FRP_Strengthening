Threshold-reliablity-based FRP strengthening of RC bridge girders
==========================

|              |                                                           |
| ------------ | --------------------------------------------------------- |
| Date:        | Sep 2015                                                  |
| Author(s):   | David Yang                                                |
| Contact:     | ynyang1988@gmail.com                                      |
| License:     | Software is released under the GNU General Public Licence.|
| Version:     | 1.0.0                                                     |

DBN-MC simulation of RC girders:
-----
* Things to be checked when switching between flexure and shear analysis:
  1. ```constants/beamConstants.py```: comment/uncomment shear dimensions for
  flexure and shear analysis, respectively
  2. ```constants/evidenceConstants.py```: appropriate CS transition years
  3. ```Life_Cycle_RC.py```: select right seed and folder to generate the samples
  and save the results

* Things to be checked when switching between analysis without and with
  evidence:
  1. ```evidence/evidence.py```: comment/uncomment CS block for analysis
  without/with evidence
  2. ```constants/evidenceConstants.py```: appropriate CS transition years for
  analysis with evidence

Reliability analysis of RC girders:
-----
* ```pyCEsmp/mainSmp.py```: change condAvailability accordingly for
  flexure, shear, independent system or correlated system
* Postprocessing functions are mostly provided in ```pyCEsmp/CEfunctions.py```,
  a small number of postprocessing functions are in ```fig/fig_source```

DBN-MC simulation of FRP-strengthened girders:
-----
* Things to be checked when switching among different strengthening schemes
  1. ```constants/beamConstants.py```: since only shear strengthening is needed,
     beam dimension for shear should be used (uncomment appropriate block)
  2. For similar reason as 1, ```constants/evidenceConstants.py``` and
     ```evidence/evidence.py``` should be checked to make sure they match with
     the analysis type: shear analysis with evidence
  3. ```constants/beamConstants.py```: comment/uncomment FRP-related dimensions
     for U-w/-anchor and U-w/o-anchor respectively
  4. ```Life_Cycle_FRP.py```: select right seed and folder to generate the samples
     and save the results
  5. ```Life_Cycle_FRP.py```: comment/uncomment lines 52 and 56 for analysis with
     and without degradation
