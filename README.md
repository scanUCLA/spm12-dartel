# spm12-dartel

Code for preprocessing of functional and structural MRI data into standardized MNI space using SPM12 and DARTEL. 

* Code in <b>/spm12_dartel_1struct</b> is for datasets that include only one structural scan (e.g. either T1 MPRAGE or T2 matched-bandwidth)
* Code in <b>/spm12_dartel_mprageMBW</b> is for datasets that include two structural scans (e.g. T1 MPRAGE *and* T2 matched-bandwidth)

<b>Instructions:</b>

Within each folder there is a <b>wrapper</b> script and a <b>run</b> function. All user-editable parameters are in an the epynomous section of the wrapper. Other sections of the wrapper script and run function shouldn't be edited unless you know what you're doing. Call only the wrapper as the wrapper will call the run function in a parfor loop. A "runStatus" struct containg each subject's pre-dartel status will be saved in the folder specified in "batchDir". The matlab workspace after pre-dartel will also be saved in "batchDir", you can use this to re-run DARTEL without re-running pre-dartel. 

<b>1struct algorithm:</b>
1) Realign functionals to mean functional (pre-dartel; parfor parallelization)
2) Co-register structural to mean functional (pre-dartel; parfor parallelization)
3) Segment structural (pre-dartel; implicitly parallaleized)
4) Create DARTEL templates
5) Normalize functionals to MNI space via DARTEL
6) Normalize structural to MNI space via DARTEL

<b>mprageMBW algorithm:</b>
1) Realign functionals to mean functional (pre-dartel; parfor parallelization)
2) Co-register MBW to mean functional (pre-dartel; parfor parallelization)
3) Co-register MPRAGE to MBW (pre-dartel; parfor parallelization)
4) Segment MPRAGE (pre-dartel; implicitly parallaleized)
5) Create DARTEL templates
6) Normalize functionals to MNI space via DARTEL
7) Normalize MPRAGE to MNI space via DARTEL
