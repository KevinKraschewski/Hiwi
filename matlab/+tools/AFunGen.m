classdef AFunGen < handle
    %AFUNGEN Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function pm = plot(this, range, varargin)
            if nargin < 2
                range = 0:.1:1000;
            elseif length(range) == 2
                range = linspace(range(1),range(2),2000);
            end
            
            [f, df] = this.getFunction;
            args = {};
            
            if ~isempty(df)
                args = {false, 1, 2};
            end
            pm = [];
            if ~isempty(varargin) && all(ishandle(varargin{1}))    
                ax = varargin{1};
                varargin(1) = [];
                for k=1:numel(ax)
                    hold(ax(k),'on');
                end
            else
                pm = PlotManager(args{:});
                pm.UseFileTypeFolders = false;
                pm.LeaveOpen = true;
                mc = metaclass(this);
                ax = pm.nextPlot(mc.Name,sprintf('Plot of %s\n%s',mc.Name,this.getConfigStr),'t [ms]','value');
                if ~isempty(df)
                    ax(2) = pm.nextPlot(mc.Name,sprintf('Plot of %s-derivative\n%s',mc.Name,this.getConfigStr),'t [ms]','value');
                end
            end
            plot(ax(1),range,f(range),varargin{:});
            if ~isempty(df)
                plot(ax(2),range,df(range),varargin{:});
            end
            if ~isempty(pm)
                pm.done;
            end
        end
        
        function plottofigure(this, fignr, range)
            if nargin < 3
                range = 0:.1:1000;
                if nargin < 2
                    fignr = gcf;
                end
            elseif length(range) == 2
                range = linspace(range(1),range(2),2000);
            end
            ax = get(fignr,'Children');
            if numel(ax) < 2
                error('Need two axes handles.');
            end
            this.plot(range, ax);
        end
        
        function sum = plus(this, other)
            % Provides an override for the simple sum of two AFunGen
            if ~isa(this,'tools.AFunGen') || ~isa(other,'tools.AFunGen')
                error('Addition not defined for non-AFunGen classes.');
            end
            sum = tools.FuncSum(this, other);
        end
    end
    
    methods(Abstract)
        [fhandle, dfhandle] = getFunction(this);
        
        str = getConfigStr(this);
    end
    
end

