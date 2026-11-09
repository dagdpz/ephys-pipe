function varargout=shadedErrorBar(x,y,errBar,lineProps,transparent)
matlabversion=datevec(version('-date'));
if nargin < 5
   transparent=0; 
end
if any(any(~isnan(errBar)))
    if matlabversion(1)>=2014
        out=shadedErrorBar_2014(x,y,errBar,'lineprops',lineProps,'transparent',transparent);
    else
        out=shadedErrorBar_2011(x,y,errBar,lineProps,transparent);
    end
else
    out=plot(x,y);
end

if ~iscell(out)
varargout={out};
else
varargout=out;
end
end
