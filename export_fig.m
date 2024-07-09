function export_fig(varargin)
matlabversion=datevec(version('-date'));
if matlabversion(1)>=2020
    exportgraphics(varargin{1}, [varargin{2} '.pdf'],'ContentType','vector') % with all respect to export_fig and previous efforts of people to make it work, let's use the existing funtionality of Matlab for saving pdfs
elseif matlabversion(1)>=2014
    export_fig_2014(varargin{:});    
else
    export_fig_2011(varargin{:});
end
end
