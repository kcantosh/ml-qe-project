function [new_data]=scale_test()

	n=[2,4,8,16,32,64,128,256,512];

	% use pca to reduce dimensionality
	% n==number pc's to retain
	% p==percentage noise	
	load anatase_rutile_runs2000.txt;
	foo=anatase_rutile_runs2000;
	[neg b]=find(foo(:,1)<-0.95);
	

	load Ox_svm2000.txt;
	load Oy_svm2000.txt;
	load Oz_svm2000.txt;
	load a_svm2000.txt;
	load c_svm2000.txt;

	load anatase_rutile_runs2000.txt
	test = anatase_rutile_runs2000;

	%[neg dum]=find(test(:,1)<0);
	neg=neg(1:500);
	f = Ox_svm2000(neg,2:513);

	Ox = Ox_svm2000(neg,1);
	Oy = Oy_svm2000(neg,1);
	Oz = Oz_svm2000(neg,1);
	aa = a_svm2000(neg,1);
	cc = c_svm2000(neg,1);

	inds=1:500;

	[a b c] = svd(f(1:500,:));
	start = 0*50+1;
	stop = (0+1)*50;
	test = start:stop;
	train = setdiff(inds,test);
	tr_y = [Oy(train), Oz(train), aa(train), cc(train)];
	test_y = [Oy(test), Oz(test), aa(test), cc(test)];
	csvwrite(['trainy_scale.txt'],tr_y);
	csvwrite(['testy_scale.txt'],test_y);

	for k=1:9
	tr_x = f(train,:)*c(:,1:n(k));
	test_x = f(test,:)*c(:,1:n(k));

	csvwrite(['trainx_',num2str(n(k)),'_scale.txt'],tr_x);
	csvwrite(['testx_',num2str(n(k)),'_scale.txt'],test_x);
	
	tic; system(["mimo-svr/bin/msvr.x trainx_",num2str(n(k)),"_scale.txt testx_",num2str(n(k)),"_scale.txt trainy_scale.txt testy_scale.txt rbf 0.8 6 0.1"]); toc;


	end
