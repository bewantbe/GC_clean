% pic typical noise distribution
% only for 2-var
% input  variables: rd,  pic_output(), fontsize, linewidth
% output variables: (none)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% residual distribution
[p, len] = size(rd);
div_num = 400;                                         % number of bins in hist
min_rd = min(min(rd,[],2));
max_rd = max(max(rd,[],2));
rd_hist_div = linspace(min_rd, max_rd, div_num);       % value range
rd_pdf = zeros(p,length(rd_hist_div));
for kk = 1:p
    rd_pdf(kk,:) = hist(rd(kk,:), rd_hist_div);
end
rd_pdf = rd_pdf / len / (max_rd - min_rd) * div_num;

%plot(rd_hist_div, rd_pdf(1,:).*rd_hist_div.^2, rd_hist_div, rd_pdf(2,:).*rd_hist_div.^2);      % 2-var

figure(12); set(gca, 'fontsize', fontsize);
%plot(rd_hist_div, log10(rd_pdf),'linewidth',linewidth);
semilogy(rd_hist_div, rd_pdf,'linewidth',linewidth);
sax=axis();  sax(3)=10^-3;  axis(sax);
xlabel('residual value');
ylabel('probability density(log10)');
hd=legend('x', 'y');
set(hd, 'fontsize',fontsize);
pic_output_color('pdf_rd_log10');
pic_data_save('pdf_rd_log10',rd_hist_div, rd_pdf);

figure(13); set(gca, 'fontsize', fontsize);
plot(rd_hist_div, rd_pdf,'linewidth',linewidth);
xlabel('residual value');
ylabel('probability density');
hd=legend('x', 'y');
set(hd, 'fontsize',fontsize);
pic_output_color('pdf_rd');
pic_data_save('pdf_rd',rd_hist_div, rd_pdf);
