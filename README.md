# tftools
Tools for analysis of time and frequency data.

MATLAB functions and scripts
----------------------------

|     |  Deviation    |
| ---- | -----|
|adev.m      |  Overlapping/non-overlapping Allan deviation, with/without gaps|
|hdev.m      | Overlapping/non-overlapping Hadamard deviation|
|mdev.m      |  Modified Allan deviation|
|tdev.m      |  Time deviation |
|theo1dev.m  |  Theo1 deviation |
|totdev.m     | Total deviation, with/without gaps|

|     | Bias correction |
| ---- | -----|
|theo1devbias.m | Bias correction for Theo1, for specified noise type|

|     | Confidence interval |
| ---- | -----|
|adevci.m  |  Estimation of confidence interval for specified noise type|
|theo1devci | ditto |

|      | Miscellaneous |
| ---- | -----|
|approxchi2inv.m | Approximate solution of the inverse chi-squared function|
|example.m | Various examples |
|freq2phase.m | Converts frequency data to phase |
|markgaps.m    | Detects and tags gaps in a time series |
|validation.m  |Runs a validation against test data from NIST|



