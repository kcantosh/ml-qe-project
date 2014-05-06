library("nnet")
ptm<-proc.time()
testInput<-read.csv('testx_2_scale.txt',header=FALSE)
trainInput<-read.csv('trainx_2_scale.txt',header=FALSE)
trainOutput<-read.csv('trainy_scale.txt',header=FALSE)
neuralNetworkModel<-nnet(trainInput,trainOutput,size=4,linout=T,skip=T,decay=5e-2,maxit=20000);
predOutput<-predict(neuralNetworkModel, testInput)
proc.time()-ptm
ptm<-proc.time()
testInput<-read.csv('testx_4_scale.txt',header=FALSE)
trainInput<-read.csv('trainx_4_scale.txt',header=FALSE)
trainOutput<-read.csv('trainy_scale.txt',header=FALSE)
neuralNetworkModel<-nnet(trainInput,trainOutput,size=4,linout=T,skip=T,decay=5e-2,maxit=20000);
predOutput<-predict(neuralNetworkModel, testInput)
proc.time()-ptm
ptm<-proc.time()
testInput<-read.csv('testx_8_scale.txt',header=FALSE)
trainInput<-read.csv('trainx_8_scale.txt',header=FALSE)
trainOutput<-read.csv('trainy_scale.txt',header=FALSE)
neuralNetworkModel<-nnet(trainInput,trainOutput,size=4,linout=T,skip=T,decay=5e-2,maxit=20000);
predOutput<-predict(neuralNetworkModel, testInput)
proc.time()-ptm
ptm<-proc.time()
testInput<-read.csv('testx_16_scale.txt',header=FALSE)
trainInput<-read.csv('trainx_16_scale.txt',header=FALSE)
trainOutput<-read.csv('trainy_scale.txt',header=FALSE)
neuralNetworkModel<-nnet(trainInput,trainOutput,size=4,linout=T,skip=T,decay=5e-2,maxit=20000);
predOutput<-predict(neuralNetworkModel, testInput)
proc.time()-ptm
ptm<-proc.time()
testInput<-read.csv('testx_32_scale.txt',header=FALSE)
trainInput<-read.csv('trainx_32_scale.txt',header=FALSE)
trainOutput<-read.csv('trainy_scale.txt',header=FALSE)
neuralNetworkModel<-nnet(trainInput,trainOutput,size=4,linout=T,skip=T,decay=5e-2,maxit=20000);
predOutput<-predict(neuralNetworkModel, testInput)
proc.time()-ptm
ptm<-proc.time()
testInput<-read.csv('testx_64_scale.txt',header=FALSE)
trainInput<-read.csv('trainx_64_scale.txt',header=FALSE)
trainOutput<-read.csv('trainy_scale.txt',header=FALSE)
neuralNetworkModel<-nnet(trainInput,trainOutput,size=4,linout=T,skip=T,decay=5e-2,maxit=20000);
predOutput<-predict(neuralNetworkModel, testInput)
proc.time()-ptm
ptm<-proc.time()
testInput<-read.csv('testx_128_scale.txt',header=FALSE)
trainInput<-read.csv('trainx_128_scale.txt',header=FALSE)
trainOutput<-read.csv('trainy_scale.txt',header=FALSE)
neuralNetworkModel<-nnet(trainInput,trainOutput,size=4,linout=T,skip=T,decay=5e-2,maxit=20000);
predOutput<-predict(neuralNetworkModel, testInput)
proc.time()-ptm
ptm<-proc.time()
testInput<-read.csv('testx_256_scale.txt',header=FALSE)
trainInput<-read.csv('trainx_256_scale.txt',header=FALSE)
trainOutput<-read.csv('trainy_scale.txt',header=FALSE)
neuralNetworkModel<-nnet(trainInput,trainOutput,size=4,linout=T,skip=T,decay=5e-2,maxit=20000);
predOutput<-predict(neuralNetworkModel, testInput)
proc.time()-ptm
ptm<-proc.time()
testInput<-read.csv('testx_512_scale.txt',header=FALSE)
trainInput<-read.csv('trainx_512_scale.txt',header=FALSE)
trainOutput<-read.csv('trainy_scale.txt',header=FALSE)
neuralNetworkModel<-nnet(trainInput,trainOutput,size=4,linout=T,skip=T,decay=5e-2,maxit=20000);
predOutput<-predict(neuralNetworkModel, testInput)
proc.time()-ptm
