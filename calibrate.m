function [out]=calibrate()

	%


	p=[0.001 0.01 0.1 ];

	load("../../trainy_0.txt");

	[ m n]=size(trainy_0);

	for k=1:3
		for i=1:30
			epsi = 0.2*i;
			for j=1:30
				C  = 0.2*j;
				
				results = zeros(1,n);
				for l=0:9
					[output text]=system(["./msvr.x ../../trainx_",num2str(l),".txt ../../testx_",num2str(l),".txt ../../trainy_",num2str(l),".txt ../../testy_",num2str(l),".txt rbf ",num2str(epsi)," ",num2str(C)," ",num2str(p(k))]);
					out((k-1)*900+ (i-1)*30+j,:,l)=[str2num(text),epsi,C,p(k)]; 
					%results += str2num(text);
				end

				%out((k-1)*900+ (i-1)*30+j,:)=[results./10,epsi,C,p(k)]; 
			end
		end
	end
