function [d]=gauge_pca(n)

% look at convergence of evalues with sample size

load Ox_svm10000.txt;
hold on;
for i=0:n

	f = Ox_svm10000(i*10+1:end,2:513);
	g=f'*f;
	%d(:,i+1) = svd(f);
	d(:,i+1) = sort(eig(g),'descend');

end	


for i=1:10

	plot(d(i,:))

end

	
