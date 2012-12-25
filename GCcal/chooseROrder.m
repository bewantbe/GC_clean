% Choose Order according to covariance
%% example:
%X = gdata(1e5, 3, 5);
%R = getcovpd(X, maxod);
%[best_od, xic] = chooseRGCOrder(R,length(X),'AIC')
%currently support 'AIC', 'BIC', 'AICc'
% if search_mode == 1, search minimal pointly

function [best_od, xic] = chooseROrder(R, len, ic_mode, search_mode)
[p, maxod] = size(R);
maxod = round(maxod/p)-1;
if ~exist('ic_mode','var')
    ic_mode = 'BIC';
end
if isempty(ic_mode)==1 || strcmpi(ic_mode,'BIC')
    fIC = @(sigma,p,m,len) log(det(sigma)) + p^2*m*log(len)/len;
elseif strcmpi(ic_mode,'AIC')
    fIC = @(sigma,p,m,len) log(det(sigma)) + 2*p^2*m/len;
elseif strcmpi(ic_mode,'AICc')
    fIC = @(sigma,p,m,len) log(det(sigma)) + 2*p^2*m/len + 2*p^2*m*(p^2*m+1)/(len-p^2*m-1)/len;
else  % BIC
    disp('ic_mode:');
    disp(ic_mode);
    error('unknown information criterion.');
end

if ~exist('search_mode','var') || search_mode == 0
    % fast mode
    minod = 3;
    mstep = 2;
    XIC = Inf(1, maxod+1);    % be careful here !
    % assume it's a convex function
    for m = minod:mstep:maxod
        [oA, sigma] = ARregression(R(:,1:m*p+p));
        XIC(m) = fIC(sigma,p,m,len);
        if (XIC(m-mstep) < XIC(m))
            break;
        end
    end
    m = m-1;
    [oA, sigma] = ARregression(R(:,1:m*p+p));
    XIC(m) = fIC(sigma,p,m,len);
    m = m-2;
    [oA, sigma] = ARregression(R(:,1:m*p+p));
    XIC(m) = fIC(sigma,p,m,len);
else
    min_ic = Inf;
    XIC = Inf(1, maxod+1);
    for m = 1:maxod
        [oA, sigma] = ARregression(R(:,1:m*p+p));
        XIC(m) = fIC(sigma,p,m,len);
        if min_ic > XIC(m)
            min_ic = XIC(m);
        end
        if XIC(m) - min_ic > 5*fIC(1,p,1,len)   % almost enough
            break;
        end
    end
end

[xic, best_od] = min(XIC);
if xic == Inf
    best_od = maxod;  % can not determine best order
end

end
