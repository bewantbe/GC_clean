id_prps = 1;
id_ps   = 1;

id_stv  = 1;
stv     = s_stv(id_stv)
slen    = round(T_segment/stv)
if (exist('signature0', 'var'))
  mlen = T_segment/stv0;
  S1 = prps_ps_stv_S  {id_prps, id_ps, id_stv}*stv*(mlen/slen/2)^2;
  fqs= prps_ps_stv_fqs{id_prps, id_ps, id_stv}*1e3;
else
  S1 = prps_ps_stv_S  {id_prps, id_ps, id_stv}*stv;
  fqs= prps_ps_stv_fqs{id_prps, id_ps, id_stv}*slen/T_segment*1e3;
end
d1 = prps_ps_stv_R(1,1, id_prps, id_ps, id_stv)
d2 = prps_ps_stv_R(2,2, id_prps, id_ps, id_stv);
dxy= prps_ps_stv_R(1,2, id_prps, id_ps, id_stv);
figure(11);
plot(fftshift(fqs), fftshift(S1(:,1,1)) - d1*stv);
%%hold on
figure(12);
plot(fftshift(fqs), fftshift(S1(:,1,2) - dxy*stv));
%%hold on
mean(S1(:, 1, 1))
max(S1(:, 1, 1))
fqs(2)-fqs(1)

%figure(11);
%plot(fftshift(fqs), fftshift(S1(:,1,1)));
%figure(12);
%plot(fftshift(fqs), fftshift(S1(:,1,2)));

%S1 = prps_ps_stv_S  {id_prps, id_ps, id_stv};
%De = prps_ps_stv_R  (1:2,1:2, id_prps, id_ps, id_stv);
%fq_cut = 250;
%[S3, fqs3] = FreqCut(S1, fqs, fq_cut);
%S3 = nuft_bias_removal(S3, De, stv, slen);
%S3 = Makeup4SpectrumFact(S3);
%figure(14);
%plot(fftshift(fqs3), fftshift(S3(:,1,1)));
%[gc, de11, de22] = getGCSapp(S3);
%gc

%id_stv = 20;
%stv     = s_stv(id_stv);
%slen    = round(T_segment/stv);
%S1 = prps_ps_stv_S  {id_prps, id_ps, id_stv}*stv;
%fqs= prps_ps_stv_fqs{id_prps, id_ps, id_stv}*slen/T_segment*1e3;
%d1 = prps_ps_stv_R  (1,1, id_prps, id_ps, id_stv);
%d2 = prps_ps_stv_R  (2,2, id_prps, id_ps, id_stv);
%dxy= prps_ps_stv_R  (1,2, id_prps, id_ps, id_stv);
%figure(11);
%hd=plot(fftshift(fqs), fftshift(S1(:,1,1))- d1*stv,'r');
%legend(hd,['stv=',num2str(stv)]);
%hold off
%figure(12);
%hd=plot(fftshift(fqs), fftshift(S1(:,1,2)-dxy*stv), 'r');
%legend(hd,['stv=',num2str(stv)]);
%hold off

%xlim([-500,500]);


%s_gc = zeros(2, length(s_stv));
%for id_stv = 1:length(s_stv);
  %stv     = s_stv(id_stv);
  %slen    = round(T_segment/stv);
  %S1 = prps_ps_stv_S  {id_prps, id_ps, id_stv};
  %fqs= prps_ps_stv_fqs{id_prps, id_ps, id_stv}*slen/T_segment*1e3;
  %De = prps_ps_stv_R  (:,:, id_prps, id_ps, id_stv);
  %S1 = nuft_bias_removal(S1, De, stv, slen);
  %[gc, de11, de22] = getGCSapp(S1);
  %s_gc(:, id_stv) = [gc(2,1), gc(1,2)];
%end

%figure(13);
%plot(s_stv, s_gc);


