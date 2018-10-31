function steadyStatePlots(par)

global Params;





x = [0.01:0.01:0.5]';


par = par2wide(par);


saveFunc = [];


for ip=1:Params.npp
    S = savingspline(par(Params.par_sind,ip));

   
    s = interp_savspline(S,x);
    saveFunc = [saveFunc s];
   

    
    
end


f1 = figure('Color',[1 1 1],'Position',[1 1 1000 700]);

plot(x,saveFunc,x,x,'LineWidth',2);

xlim([-0.02 0.5]);

