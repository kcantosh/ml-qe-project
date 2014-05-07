ml-qe-project
=============

dft &amp; machine learning for materials characterization

#step 1; generate some data
the complete workflow for doing this is in here : batch_anatase2000.sh, as run on stampede (TACC). YMMV, but you will need at least the pseudo potentials in this directory for Ti &O, reconstructed for GIPAW, and a compilation of the MAS simulation code in mqmas*c, compiled as sim_mas.x, and of course pw.x/gipaw.x as found in quantum espresso. The output is in *svm* files, one per unit cell parameter.

#step 2; scaling tests
after generating data, or using the data available in this repos, you can test the scaling of MSVR or NNET in R. To do the former, simply run scale_test.m in octave, after cloning and building https://github.com/wjb19/mimo-svr.git in this directory. Similarly, you can test the ANN after running script in octave (which generates test sets of dimensions 2,4,8,16,32,64,128,256,512) by running:
R --no-save < nnet_scale.R

#step 3; pretty/process data
pca_comp.m will allow you to compress the data and create 10 data folds for training/testing machine learning methods. Check out the file for more details but basically all you need to do is decide how many principal components you want to retain, and what threshold for Cq. Note that the latter will have impact on the final dataset you select ie., agressive thresholding could easily leave you with much less than 500 examples to work with, implicit to this script.

#step 4; machine learning
msvr (multiple-input, multiple-output SVR) will need to be calibrated, particularly the constants epsilon and C need tuning. This is all accomplished using the calibrate.m script eg., run in octave thusly:
out =calibrate;

For the ANN :
R --no-save < nnet_script.R
then cleanup the files for octave:
./clean_comp.sh
and check results in octave with nnet_validate; should be four columns of output, one row per datafold (10 total). 
