
load pred0_clean.txt;
load pred1_clean.txt;
load pred2_clean.txt;
load pred3_clean.txt;
load pred4_clean.txt;
load pred5_clean.txt;
load pred6_clean.txt;
load pred7_clean.txt;
load pred8_clean.txt;
load pred9_clean.txt;

load testy_0.txt;
load testy_1.txt;
load testy_2.txt;
load testy_3.txt;
load testy_4.txt;
load testy_5.txt;
load testy_6.txt;
load testy_7.txt;
load testy_8.txt;
load testy_9.txt;

ln = length(pred0_clean);

y(1,:) = sqrt(sum((pred0_clean-testy_0).^2) / ln);
y(2,:) = sqrt(sum((pred1_clean-testy_1).^2) / ln);
y(3,:) = sqrt(sum((pred2_clean-testy_2).^2) / ln);
y(4,:) = sqrt(sum((pred3_clean-testy_3).^2) / ln);
y(5,:) = sqrt(sum((pred4_clean-testy_4).^2) / ln);
y(6,:) = sqrt(sum((pred5_clean-testy_5).^2) / ln);
y(7,:) = sqrt(sum((pred6_clean-testy_6).^2) / ln);
y(8,:) = sqrt(sum((pred7_clean-testy_7).^2) / ln);
y(9,:) = sqrt(sum((pred8_clean-testy_8).^2) / ln);
y(10,:) = sqrt(sum((pred9_clean-testy_9).^2) / ln);
