This file contains the data used to fit the HMMs (regularized at 12-hour resolution using CTMM and split into bursts any time there was a gap >24 hours). Columns include:

ID: this is the identifier of the bursts used to fit the HMM. It's a combination of "deploy_ID" and "burst"
uniqueID: animal identifier
collarID: collar identifier
...
[self-explanatory columns]
...
step: step length calculated by momentuHMM via prepData() function
angle: turn angle calculated by momentuHMM via prepData() function
deploy_ID: combination of animal and collar ID
t_: timestamp (Mountain Time)
uncertainty: whether the location was recorded ("None") or interpolated ("Uncertainty")
data_gap: whether we had any data at all on this date
burst_: split each time there is a data gap
x: UTM Easting
y: UTM Northing
vit_2states: Viterbi assignments from the 2-state HMM
vit_3states: Viterbi assignments from the 3-state HMM