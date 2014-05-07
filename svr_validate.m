load msvr_calibration.mat

for i=1:10
tot=out(:,1,i)+out(:,2,i)+out(:,3,i)+out(:,4,i)+out(:,5,i);

[a b]=min(tot);

y(i,:) = out(b,:,i);

%[output text]=system(["mimo-svr/bin/msvr.x trainx_",num2str(i-1),".txt testx_",num2str(i-1),".txt trainy_",num2str(i-1),".txt testy_",num2str(i-1),".txt rbf ",num2str(y(i,5))," ",num2str(y(i,6))," ",num2str(y(i,7))]);

%results(i,:)=str2num(text);

end

