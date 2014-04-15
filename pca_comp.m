function [new_data]=pca_comp(n)

	% use pca to reduce dimensionality
	% n==number pc's to retain
	% p==percentage noise	
	load anatase_rutile_runs2000.txt;
	foo=anatase_rutile_runs2000;
	[neg b]=find(foo(:,1)<-1.5);
	

	load Ox_svm2000.txt;
	load Oy_svm2000.txt;
	load Oz_svm2000.txt;
	load a_svm2000.txt;
	load c_svm2000.txt;

	load anatase_rutile_runs2000.txt
	test = anatase_rutile_runs2000;

	%[neg dum]=find(test(:,1)<0);
	neg=neg(1:300);
	f = Ox_svm2000(neg,2:513);

	Ox = Ox_svm2000(neg,1);
	Oy = Oy_svm2000(neg,1);
	Oz = Oz_svm2000(neg,1);
	aa = a_svm2000(neg,1);
	cc = c_svm2000(neg,1);

	inds=1:300;

	[a b c] = svd(f(1:300,:));

	for k=0:9

	start = k*30+1;
	stop = (k+1)*30;
	test = start:stop;
	train = setdiff(inds,test);


	%retain n pc's
	tr_x = f(train,:)*c(:,1:n);
	test_x = f(test,:)*c(:,1:n);

	tr_y = [Ox(train), Oy(train), Oz(train), aa(train), cc(train)];
	test_y = [Ox(test), Oy(test), Oz(test), aa(test), cc(test)];


	csvwrite(['trainx_',num2str(k),'.txt'],tr_x);
	csvwrite(['testx_',num2str(k),'.txt'],test_x);
	csvwrite(['trainy_',num2str(k),'.txt'],tr_y);
	csvwrite(['testy_',num2str(k),'.txt'],test_y);


	


%for i=1:length(train)
	%	for j=1:n-1
	%	tmp=[num2str(tr_x(i,j)),','];
	%	end
	%	tmp=[tmp,num2str(tr_x(i,n))];
	%
	%	fprintf(fid1,'%s \n',tmp);

	%	end

	%	fclose(fid1);
%for i=1:length(test)
	%	for j=1:n-1
	%	tmp=[num2str(test_x(i,j)),','];
	%	end
	%	tmp=[tmp,num2str(test_x(i,n))];
	%
	%	fprintf(fid2,'%s \n',tmp);
	%
	%	end

	%	fclose(fid2);





	end
